#!/usr/bin/env python3
import argparse
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import List, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
ELF_PATH = SCRIPT_DIR / "obj" / "mayo_P3_fault.elf"
RESULTS_DIR = SCRIPT_DIR / "bash_script_results"
RAW_LOG = RESULTS_DIR / "instruction_skip_raw.txt"
HITS_FILE = RESULTS_DIR / "instruction_skip_hits.txt"


class Instruction:
    def __init__(self, address: int, bytes_: List[int], asm: str):
        self.address = address
        self.bytes = bytes_
        self.asm = asm

    def __repr__(self) -> str:
        return f"0x{self.address:x}: {self.bytes} {self.asm}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Find instruction-skip faults for an ARM ELF under QEMU/GDB")
    parser.add_argument(
        "--symbol",
        default="P1_times_O",
        help="Function name or substring to disassemble (default: P1_times_O)",
    )
    parser.add_argument(
        "--address",
        default=None,
        help="Optional single instruction address to test, e.g. 0x3324",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=25,
        help="Per-instruction GDB timeout in seconds",
    )
    return parser.parse_args()


def find_symbol(elf_path: Path, symbol_substring: str) -> str:
    result = subprocess.run(
        ["nm", "-n", str(elf_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"nm failed: {result.stdout}")

    matches: List[str] = []
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    for line in lines:
        parts = line.split()
        if len(parts) < 3:
            continue
        sym = parts[-1]
        if symbol_substring in sym:
            matches.append(sym)

    if not matches:
        raise RuntimeError(f"Could not find symbol matching {symbol_substring!r}")

    exact_matches = [sym for sym in matches if sym == symbol_substring]
    if exact_matches:
        return exact_matches[0]

    non_constprop_matches = [sym for sym in matches if ".constprop." not in sym]
    if non_constprop_matches:
        return non_constprop_matches[0]

    return matches[0]


def get_disassembly(elf_path: Path, symbol_name: str) -> List[Instruction]:
    candidates = [
        ["arm-none-eabi-objdump", "-d", "-w", "-C", f"--disassemble={symbol_name}", str(elf_path)],
        ["objdump", "-d", "-w", "-C", f"--disassemble={symbol_name}", str(elf_path)],
        ["objdump", "-d", "-w", "-m", "arm", f"--disassemble={symbol_name}", str(elf_path)],
    ]

    last_error = None
    for cmd in candidates:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if result.returncode == 0 and result.stdout.strip():
            output = result.stdout
            break
        last_error = result.stdout
    else:
        raise RuntimeError(f"objdump failed: {last_error}")

    output = result.stdout

    instructions: List[Instruction] = []
    for raw_line in output.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("Disassembly"):
            continue

        if line.startswith("<") and line.endswith(":"):
            continue

        if not re.search(r"^[0-9a-fA-F]+:\s", line):
            continue

        try:
            prefix, rest = line.split(":", 1)
            addr = int(prefix.strip(), 16)
        except ValueError:
            continue

        remainder = rest.strip()
        if not remainder:
            continue

        parts = remainder.split("\t")
        if not parts:
            continue

        byte_field = parts[0].strip()
        byte_tokens = byte_field.split()
        instruction_bytes: List[int] = []
        valid_bytes = True

        for token in byte_tokens:
            if len(token) % 2 != 0:
                valid_bytes = False
                break
            try:
                if len(token) == 2:
                    instruction_bytes.append(int(token, 16))
                elif len(token) == 4:
                    instruction_bytes.extend([
                        int(token[0:2], 16),
                        int(token[2:4], 16),
                    ])
                else:
                    for idx in range(0, len(token), 2):
                        instruction_bytes.append(int(token[idx:idx + 2], 16))
            except ValueError:
                valid_bytes = False
                break

        if valid_bytes and instruction_bytes:
            asm = " ".join(part.strip() for part in parts[1:] if part.strip())
        else:
            asm = remainder

        instructions.append(Instruction(addr, instruction_bytes, asm))

    return instructions


def build_gdb_script(elf_path: Path, instruction: Instruction, runtime_addr: int, skip_bytes: int) -> str:
    def thumb_addr(addr: int) -> int:
        return addr | 1

    script = []
    script.append("set pagination off")
    script.append("set confirm off")
    script.append("set architecture arm")
    script.append(f"file {elf_path}")
    script.append("target remote :1234")
    script.append(f"break *0x{thumb_addr(runtime_addr):x}")
    script.append("continue")
    script.append(f"set $pc = 0x{thumb_addr(runtime_addr + skip_bytes):x}")
    script.append("continue")
    script.append("quit")
    return "\n".join(script) + "\n"


def run_qemu_and_gdb(elf_path: Path, instruction: Instruction, timeout_s: int) -> Tuple[str, bool]:
    qemu_cmd = [
        "qemu-system-arm",
        "-M",
        "mps2-an386",
        "-kernel",
        str(elf_path),
        "-nographic",
        "-semihosting",
        "-S",
        "-gdb",
        "tcp::1234",
    ]

    qemu_proc = subprocess.Popen(
        qemu_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    try:
        gdb_output = ""
        gdb_ok = False
        deadline = time.time() + 10
        # For this bare-metal ELF, the addresses from objdump/nm are already the
        # runtime VMA addresses used by QEMU/GDB. Do not add an arbitrary base.
        runtime_addr = instruction.address
        while time.time() < deadline:
            # Use the actual instruction width parsed from objdump. Thumb-2
            # instructions such as BL can be 4 bytes wide, so skipping by 2 bytes
            # would leave execution in the middle of the instruction.
            skip_length = len(instruction.bytes) if instruction.bytes else 2
            gdb_script = build_gdb_script(elf_path, instruction, runtime_addr, skip_length)
            with tempfile.NamedTemporaryFile("w", delete=False) as handle:
                handle.write(gdb_script)
                script_path = handle.name

            try:
                proc = subprocess.run(
                    ["gdb-multiarch", "-q", "-batch", "-x", script_path],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=timeout_s,
                    check=False,
                )
                gdb_output = proc.stdout
                gdb_ok = proc.returncode == 0
            except subprocess.TimeoutExpired:
                gdb_output = "[gdb timed out]"
                gdb_ok = False
            finally:
                try:
                    os.remove(script_path)
                except FileNotFoundError:
                    pass

            if gdb_ok or "Remote debugging using" in gdb_output or "Cannot connect" not in gdb_output:
                break
            time.sleep(0.5)

        qemu_stdout = ""
        try:
            qemu_stdout, _ = qemu_proc.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            qemu_proc.kill()
            qemu_stdout, _ = qemu_proc.communicate(timeout=5)

        combined = (gdb_output + "\n" + qemu_stdout).strip()
        passed = "PASS" in combined
        return combined, passed
    finally:
        if qemu_proc.poll() is None:
            qemu_proc.kill()
            qemu_proc.wait(timeout=5)


def write_results(results: List[Tuple[Instruction, str]], hits_file: Path, raw_log: Path) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with raw_log.open("w", encoding="utf-8") as handle:
        for instruction, output in results:
            handle.write(f"========================================\n")
            handle.write(f"0x{instruction.address:x}  {instruction.asm}\n")
            handle.write(f"orig_bytes={instruction.bytes}\n")
            handle.write(output)
            handle.write("\n")

    with hits_file.open("a", encoding="utf-8") as handle:
        for instruction, output in results:
            if "PASS" in output:
                handle.write(f"0x{instruction.address:x}  {instruction.asm}\n")


def main() -> None:
    args = parse_args()

    if not ELF_PATH.exists():
        raise SystemExit(f"ELF not found: {ELF_PATH}")

    symbol = find_symbol(ELF_PATH, args.symbol)
    instructions = get_disassembly(ELF_PATH, symbol)
    if not instructions:
        raise SystemExit("No instructions were parsed from the target function")

    if args.address is not None:
        target_addr = int(args.address, 16)
        instructions = [insn for insn in instructions if insn.address == target_addr]
        if not instructions:
            raise SystemExit(f"No instruction found at address {args.address}")

    print(f"Disassembled {len(instructions)} instructions from {symbol}")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    raw_results: List[Tuple[Instruction, str]] = []

    for idx, instruction in enumerate(instructions, 1):
        print(f"[{idx}/{len(instructions)}] Trying skip at 0x{instruction.address:x}: {instruction.asm}")
        output, passed = run_qemu_and_gdb(ELF_PATH, instruction, args.timeout)
        raw_results.append((instruction, output))
        if passed:
            print(f"[HIT] 0x{instruction.address:x}: PASS detected")

    write_results(raw_results, HITS_FILE, RAW_LOG)
    print(f"Raw log: {RAW_LOG}")
    print(f"Hits file: {HITS_FILE}")


if __name__ == "__main__":
    main()
