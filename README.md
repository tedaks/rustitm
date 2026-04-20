# RustITM — Rust Port of the ITS Irregular Terrain Model

A Rust implementation of the ITS Irregular Terrain Model (ITM) for radio wave propagation prediction. ITM predicts terrestrial radio wave attenuation for frequencies between 20 MHz and 20 GHz, based on electromagnetic theory and empirical models by Anita Longley and Phil Rice. Propagation mechanisms include free space loss, diffraction, and troposcatter.

This crate is a direct port of the [NTIA/ITM C++ library](https://github.com/NTIA/itm), version 1.3 (functionally identical to FORTRAN version 1.2.2).

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
rustitm = "0.1"
```

### Example — Area Mode

```rust
use rustitm::entry::itm_area_tls;

let result = itm_area_tls(
    15.0, 3.0,           // h_tx, h_rx (meters)
    2, 1,                // climate, polarization
    100.0,               // d__km
    5.0,                 // delta_h__meter
    4,                   // tx_siting_criteria
    301.0,               // N_0
    980.0,               // f__mhz
    1,                   // pol
    15.0,                // epsilon
    0.008,               // sigma
    2,                   // mdvar
    92.0, 53.0, 97.0,    // time, location, situation
).unwrap();

println!("Basic transmission loss: {:.2} dB", result.a__db);
```

### Example — Point-to-Point Mode

```rust
use rustitm::entry::itm_point_to_point;

let result = itm_point_to_point(
    100.0, 50.0,         // h_tx, h_rx (meters)
    3, 1,                // climate, polarization
    301.0,               // N_0
    1000.0,              // f__mhz
    1,                   // pol
    15.0,                // epsilon
    0.008,               // sigma
    3,                   // mdvar
    &pfl,                // terrain profile
    50.0, 80.0, 43.0,    // time, location, situation
).unwrap();
```

## Inputs

ITM operates in two prediction modes: **Area** and **Point-to-Point**. Variability can be specified via either time/location/situation or confidence/reliability.

### Common Inputs

| Variable          | Type     | Units | Limits                  | Description                          |
|-------------------|----------|-------|-------------------------|--------------------------------------|
| `h_tx__meter`     | `f64`    | m     | 0.5 – 3000              | Structural height of TX              |
| `h_rx__meter`     | `f64`    | m     | 0.5 – 3000              | Structural height of RX              |
| `climate`        | `i32`    |       | 1–7                     | Radio climate (see below)           |
| `polarization`    | `i32`    |       | 0, 1                    | 0 = Horizontal, 1 = Vertical        |
| `N_0`             | `f64`    | N-Units | 250 – 400             | Minimum monthly mean surface refractivity (sea level) |
| `f__mhz`          | `f64`    | MHz   | 20 – 20000              | Frequency                            |
| `epsilon`         | `f64`    |       | > 1                     | Relative permittivity                |
| `sigma`           | `f64`    | S/m   | > 0                     | Conductivity                         |
| `mdvar`           | `i32`    |       | 0–3 (+10, +20 options)  | Mode of variability                  |

**Climate codes:** 1 = Equatorial, 2 = Continental Subtropical, 3 = Maritime Subtropical, 4 = Desert, 5 = Continental Temperate, 6 = Maritime Temperate Over Land, 7 = Maritime Temperate Over Sea.

**Variability modes:** 0 = Single Message, 1 = Accidental, 2 = Mobile, 3 = Broadcast. Add +10 to eliminate location variability, +20 to eliminate situation variability.

### Area Mode Additional Inputs

| Variable             | Type | Units | Limits     | Description                      |
|----------------------|------|-------|------------|----------------------------------|
| `d__km`              | `f64` | km   | > 0        | Path distance                    |
| `delta_h__meter`     | `f64` | m    | ≥ 0        | Terrain irregularity parameter   |
| `tx_siting_criteria` | `i32` |       | 0–2        | TX siting: 0=Random, 1=Careful, 2=Very Careful |
| `rx_siting_criteria` | `i32` |       | 0–2        | RX siting: 0=Random, 1=Careful, 2=Very Careful |

### Point-to-Point Mode Additional Inputs

| Variable | Type        | Units | Description                                        |
|----------|-------------|-------|----------------------------------------------------|
| `pfl`    | `&[f64]`    | m     | Terrain profile in PFL format: `[n, res, h0, h1, ...]` where `n` is point count − 1 and `res` is resolution in meters |

### Variability Inputs

Provide either time/location/situation values **or** confidence/reliability values:

| Variable     | Type   | Limits          | Description              |
|--------------|--------|-----------------|--------------------------|
| `time`       | `f64`  | 0 – 100         | Time variability         |
| `location`   | `f64`  | 0 – 100         | Location variability     |
| `situation` | `f64`  | 0 – 100         | Situation variability    |
| `confidence` | `f64`  | 0 – 100         | Confidence               |
| `reliability`| `f64`  | 0 – 100         | Reliability              |

## Outputs

| Variable     | Type     | Units | Description                    |
|--------------|----------|-------|--------------------------------|
| `a__db`      | `f64`    | dB     | Basic transmission loss        |
| `warnings`   | `i32`    |        | Warning flags (bitmask)       |
| `inter_values` | `IntermediateValues` | | Intermediate values (via `_Ex` variants) |

## Intermediate Values

Functions suffixed with `_Ex` return a struct containing intermediate calculation values:

| Variable          | Type      | Units     | Description                              |
|-------------------|-----------|-----------|------------------------------------------|
| `theta_hzn`       | `[f64; 2]` | radians  | Terminal horizon angles                  |
| `d_hzn__meter`    | `[f64; 2]` | m        | Terminal horizon distances              |
| `h_e__meter`      | `[f64; 2]` | m        | Effective terminal heights               |
| `n_s`             | `f64`     | N-Units   | Surface refractivity                    |
| `delta_h__meter`  | `f64`     | m         | Terrain irregularity parameter          |
| `a_ref__db`       | `f64`     | dB        | Reference attenuation                   |
| `a_fs__db`        | `f64`     | dB        | Free space basic transmission loss      |
| `d__km`           | `f64`     | km        | Path distance                           |
| `mode`            | `i32`     |           | Propagation mode: 1=LoS, 2=Diffraction, 3=Troposcatter |

## Error Handling

Functions return `Result<T, ItmError>`. A complete list of error and warning codes is documented in the source.

## Validation

The crate is validated against the NTIA/ITM reference vectors. All tests pass within 0.1 dB tolerance:

```
$ cargo test
   Compiling rustitm v0.1.0
    Finished test [unoptimized + debuginfo] target(s) in 0.49s
     Running tests/reference.rs
test area_tls_csv_matches_ref ... ok
test p2p_tls_cmd_example_matches_ref ... ok

test result: ok. 2 passed; 0 failed; 0 ignored
```

**Area mode** — 5 rows from `area.csv` (NTIA/ITM reference vectors): all match within 0.1 dB.

**Point-to-Point mode** — `cmd_examples` from the NTIA/ITM C++ repository (`i_p2ptls.txt` + `pfl.txt` → expected 114.5 dB): matches within 0.1 dB.

## Configure and Build

```sh
cargo build --release
```

The crate compiles with Rust 1.56+ (edition 2021). LTO is enabled for release builds.

## References

* G.A. Hufford, A.G. Longley, W.A. Kissick, [_A Guide to the Use of the ITS Irregular Terrain Model in the Area Prediction Mode_](https://www.its.bldrdoc.gov/publications/details.aspx?pub=2091), NTIA Technical Report TR-82-100, April 1982.
* G.A. Hufford, [_The ITS Irregular Terrain Model, version 1.2.2 Algorithm_](https://www.its.bldrdoc.gov/media/50676/itm_alg.pdf).
* G.A. Hufford, [_1985 ITM Memo_](https://www.its.bldrdoc.gov/media/50675/Hufford_1985_Memo.pdf).
* G.A. Hufford, [_The Irregular Terrain Model_](https://www.its.bldrdoc.gov/media/50674/itm.pdf).
* A.G. Longley and P.L. Rice, [_Prediction of Tropospheric Radio Transmission Loss Over Irregular Terrain: A Computer Method - 1968_](https://www.its.bldrdoc.gov/publications/details.aspx?pub=2784), NTIA Technical Report ERL 79-ITS 67, July 1968.

## License

MIT
