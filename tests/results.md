# Validation Results

Validated against the NTIA/ITS C++ ITM reference in `/home/bortre/02-lab/sources/itm`.

- Reference C++ built from `itm/src` + `itm/include` on Linux (g++ -O2 -std=c++17).
- Rust tests: `cargo test --release` (see `tests/reference.rs`).
- Diagnostic: `cargo run --release --example diag` (see `examples/diag.rs`).

## Final status: PASS (after fix)

```
running 2 tests
test p2p_tls_cmd_example_matches_ref ... ok
test area_tls_csv_matches_ref ... ok

test result: ok. 2 passed; 0 failed
```

## p2p_tls — cmd_examples (itm/cmd_examples/i_p2ptls.txt + pfl.txt)

Reference (o_p2ptls.txt): **114.5 dB**
Rust: within 0.1 dB. PASS.

## area_tls — itm/area.csv

Tolerance: 0.1 dB (reference values reported to 1 decimal).

| Row | h_tx | h_rx | dh   | mdvar | d_km | tx | rx | f_mhz | pol | climate | t/loc/sit | expected | rust     | cpp_ref  |
|-----|------|------|------|-------|------|----|----|-------|-----|---------|-----------|----------|----------|----------|
| 2   | 10   | 1    | 0    | 0     | 16   | 0  | 0  | 230   | 0   | 5       | 87/50/50  | 152.5    | 152.5052 | 152.5052 |
| 3   | 3    | 1.5  | 10   | 1     | 10   | 1  | 0  | 450   | 0   | 5       | 40/28/25  | 133.0    | 133.0482 | 133.0482 |
| 4   | 15   | 3    | 5    | 2     | 100  | 2  | 1  | 980   | 1   | 4       | 92/53/97  | 224.1    | 224.0971 | 224.0971 |
| 5   | 3    | 5    | 20   | 3     | 75   | 0  | 1  | 3100  | 1   | 2       | 50/80/43  | 205.1    | 205.0695 | 205.0695 |
| 6   | 1.5  | 10   | 45   | 0     | 25   | 1  | 2  | 8900  | 1   | 1       | 70/99/26  | 156.0    | 155.9746 | 155.9746 |

All 5 rows agree with both the reference CSV and the C++ reference implementation to within 0.05 dB.

## Pre-fix state (for reference)

Rows 4 and 5 previously diverged:

| Row | expected | rust (pre-fix) | diff   |
|-----|----------|----------------|--------|
| 4   | 224.1    | 222.4641       | -1.64  |
| 5   | 205.1    | 204.1262       | -0.97  |

Intermediate-value dump confirmed the gap was entirely in `A_ref__db` with
`mode=3` (troposcatter) on both failing rows; all other intermediates
(`h_e`, `d_hzn`, `theta_hzn`, `N_s`, `delta_h`, `A_fs`) matched C++ exactly.

### Root cause

`src/troposcatter.rs`, inside the `eta_s < 1.0` interpolation branch.

C++ (TroposcatterLoss.cpp:103):
```cpp
H_0 = eta_s*H_0 + (1-eta_s)*10*log10(pow((1+SQRT2/r_1)*(1+SQRT2/r_2),2)
                                     * (r_1+r_2)/(r_1+r_2+2*SQRT2));
```

Rust (pre-fix):
```rust
h_0_result = eta_s * h_0_result + (1.0-eta_s) * 10.0
    * ((1.0+SQRT2/r_1)*(1.0+SQRT2/r_2)).powi(2)
    * (r_1+r_2) / (r_1+r_2+2.0*SQRT2).log10();
```

Rust method-call precedence means `x.log10()` binds tightly, so this
parsed as `... * (r_1+r_2) / log10(r_1+r_2+2*SQRT2)` — the `log10` wrapped
only the denominator instead of the full product.

### Fix

```rust
h_0_result = eta_s * h_0_result + (1.0 - eta_s) * 10.0
    * (((1.0 + SQRT2 / r_1) * (1.0 + SQRT2 / r_2)).powi(2)
       * (r_1 + r_2) / (r_1 + r_2 + 2.0 * SQRT2)).log10();
```

## Coverage

Tested: `itm_p2p_tls`, `itm_area_tls` (TLS variants).

Not covered by these tests:
- `itm_p2p_cr`, `itm_area_cr` (CR variants) — same underlying code path,
  TLS coverage is indicative.
- Error-return and warning-flag paths.
- `itm/p2p.csv` — 5 rows, but the CSV does not couple cases to specific
  entries in `itm/pfls.csv`, so no unambiguous per-row expected value.
