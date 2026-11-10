#!/usr/bin/env python3
"""Generate deterministic LP instances that isolate presolve reductions.

The files are intentionally exact in binary floating point (small integer
coefficients and scale factors 2, 3, -1) so the parallel-row rule can require
bit-exact proportionality without accepting approximate duplicates.
"""

from pathlib import Path
import random

ROOT = Path(__file__).resolve().parent / "data"
RNG = random.Random(7)


def fmt(value: float) -> str:
    return f"{value:.17g}"


def write_lp(name: str, sense: str, objective: list[float], rows: list[tuple],
             bounds: list[tuple], generals: list[int] | None = None,
             binaries: list[int] | None = None) -> None:
    lines = [sense, " obj:" + "".join(
        f" {fmt(c)} x{j + 1}" for j, c in enumerate(objective) if c != 0.0
    ), "Subject To"]
    for index, (lower, upper, terms) in enumerate(rows, start=1):
        expression = "".join(f" {fmt(a)} x{j + 1}" for j, a in terms if a != 0.0)
        if lower == upper:
            lines.append(f" r{index}:{expression} = {fmt(lower)}")
        elif lower is None:
            lines.append(f" r{index}:{expression} <= {fmt(upper)}")
        elif upper is None:
            lines.append(f" r{index}:{expression} >= {fmt(lower)}")
        else:
            lines.append(f" r{index}: {fmt(lower)} <={expression} <= {fmt(upper)}")
    if bounds:
        lines.append("Bounds")
        for j, (lower, upper) in enumerate(bounds, start=1):
            if lower == upper:
                lines.append(f" x{j} = {fmt(lower)}")
            else:
                lower_text = "-inf" if lower is None else fmt(lower)
                upper_text = "+inf" if upper is None else fmt(upper)
                lines.append(f" {lower_text} <= x{j} <= {upper_text}")
    if binaries:
        lines.append("Binary")
        for start in range(0, len(binaries), 20):
            lines.append(" " + " ".join(f"x{j + 1}" for j in binaries[start:start + 20]))
    if generals:
        lines.append("Generals")
        for start in range(0, len(generals), 20):
            lines.append(" " + " ".join(f"x{j + 1}" for j in generals[start:start + 20]))
    lines.append("End")
    (ROOT / name).write_text("\n".join(lines) + "\n", encoding="ascii")


def forcing_lp() -> None:
    """150 three-column forcing chains plus a 300-variable core LP."""
    objective: list[float] = []
    bounds: list[tuple] = []
    rows: list[tuple] = []
    for block in range(20):
        a, b, c = (block * 3 + offset for offset in range(3))
        # a + b = 20 with both in [0,10] fixes a = b = 10.
        rows.append((20.0, 20.0, [(a, 1.0), (b, 1.0)]))
        # b + c >= 15 with b = 10 and c <= 5 fixes c = 5.
        rows.append((15.0, None, [(b, 1.0), (c, 1.0)]))
        objective.extend([1.0, 2.0, 3.0])
        bounds.extend([(0.0, 10.0), (0.0, 10.0), (0.0, 5.0)])

    core_start = len(objective)
    core_cols = 40
    core_rows = 20
    for _ in range(core_cols):
        objective.append(float(-RNG.randint(1, 9)))
        bounds.append((0.0, 10.0))
    for row in range(core_rows):
        a = core_start + 2 * row
        b = a + 1
        # x + y <= 21 leaves x = y = 10 feasible without forcing equality.
        rows.append((None, 21.0, [(a, 1.0), (b, 1.0)]))
    write_lp("forcing_lp.lp", "Minimize", objective, rows, bounds)


def parallel_lp() -> None:
    """200 base rows, each repeated with exact scale factors 2, 3, and -1."""
    columns = 20
    objective = [float(-RNG.randint(1, 9)) for _ in range(columns)]
    bounds = [(0.0, 10.0) for _ in range(columns)]
    rows: list[tuple] = []
    for _ in range(10):
        chosen = RNG.sample(range(columns), 2)
        base = [(column, float(RNG.randint(1, 9))) for column in chosen]
        upper = float(sum(coefficient * 10.0 for _, coefficient in base) * 0.8)
        rows.append((None, upper, base))
        rows.append((None, 2.0 * upper, [(column, 2.0 * coefficient) for column, coefficient in base]))
        rows.append((None, 3.0 * upper, [(column, 3.0 * coefficient) for column, coefficient in base]))
        # -1 * activity <= -upper  <=>  activity >= upper
        rows.append((-upper, None, [(column, -coefficient) for column, coefficient in base]))
    write_lp("parallel_lp.lp", "Minimize", objective, rows, bounds)


def mip_forcing() -> None:
    """50 forcing binary blocks and 50 random binary knapsacks."""
    columns = 120
    objective = [float(RNG.randint(1, 9)) for _ in range(columns)]
    bounds = [(0.0, 1.0) for _ in range(columns)]
    rows: list[tuple] = []
    for block in range(20):
        chosen = list(range(block * 4, block * 4 + 4))
        rows.append((4.0, 4.0, [(column, 1.0) for column in chosen]))
    for _ in range(20):
        chosen = RNG.sample(range(80, columns), 10)
        terms = [(column, float(RNG.randint(1, 9))) for column in chosen]
        rows.append((None, float(sum(coefficient for _, coefficient in terms) * 0.7), terms))
    write_lp("mip_forcing.lp", "Minimize", objective, rows, bounds,
             binaries=list(range(columns)))


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    forcing_lp()
    parallel_lp()
    mip_forcing()
    for path in sorted(ROOT.iterdir()):
        print(f"{path.name}: {path.stat().st_size} bytes")


if __name__ == "__main__":
    main()
