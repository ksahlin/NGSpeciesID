//! Locating the reference interpreter, for the oracles that replay CPython.
//!
//! # Three bugs this file exists to fix, all carried across from the isONclust port
//!
//! `phred_oracle`, `pyfloat_oracle` and `pyround_oracle` each had their own copy
//! of this logic, and every copy:
//!
//! 1. **defaulted to `~/miniforge3/envs/isonclust-ref`** — the *other* project's
//!    environment. It exists on the machine this port was written on, so the
//!    three oracles ran, passed, and were validating against an interpreter
//!    that is not this repository's reference. They were green for the wrong
//!    reason.
//! 2. **panicked when it was absent**, so `cargo test` failed outright on any
//!    clean clone. CI found that on all four targets; a contributor would have
//!    found it on their first checkout.
//! 3. named its escape hatch `ISONCLUST_SKIP_PYTHON_ORACLE`.
//!
//! # The contract now
//!
//! * `REF_PYTHON` wins, as it does everywhere else in this repository.
//! * Otherwise `~/miniforge3/envs/ngspeciesid-ref/bin/python`, which is what
//!   `bench/setup_reference_env.sh` builds.
//! * Absent → **skip, loudly**, naming what went unverified. A test suite that
//!   cannot run on a clean clone is a test suite people stop running.
//! * `NGSPECIESID_SKIP_PYTHON_ORACLE=1` skips explicitly.
//! * `NGSPECIESID_REQUIRE_PYTHON_ORACLE=1` turns a skip back into a failure.
//!   **CI sets this** in the job that has the environment — otherwise "skip"
//!   quietly becomes "never runs anywhere", which is the failure mode this
//!   whole file is about.

/// The reference interpreter, or `None` if these oracles must be skipped here.
pub fn reference_python() -> Option<String> {
    let required = std::env::var("NGSPECIESID_REQUIRE_PYTHON_ORACLE").is_ok();

    if std::env::var("NGSPECIESID_SKIP_PYTHON_ORACLE").is_ok() {
        assert!(
            !required,
            "NGSPECIESID_SKIP_PYTHON_ORACLE and NGSPECIESID_REQUIRE_PYTHON_ORACLE are both set"
        );
        eprintln!("  SKIPPED by NGSPECIESID_SKIP_PYTHON_ORACLE -- this oracle is unverified here");
        return None;
    }

    let py = std::env::var("REF_PYTHON").unwrap_or_else(|_| {
        format!(
            "{}/miniforge3/envs/ngspeciesid-ref/bin/python",
            std::env::var("HOME").unwrap_or_default()
        )
    });

    if std::path::Path::new(&py).is_file() {
        return Some(py);
    }

    let msg = format!(
        "the reference interpreter is not at {py}.\n\
         Build it with bench/setup_reference_env.sh, or point REF_PYTHON at one."
    );
    assert!(
        !required,
        "NGSPECIESID_REQUIRE_PYTHON_ORACLE is set but {msg}"
    );
    eprintln!("  SKIPPED -- {msg}");
    None
}
