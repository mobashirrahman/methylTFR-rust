use anyhow::Result;
use clap::Parser;

use methyltfr::cli::{Cli, main_with};

/// anyhow is used here and nowhere else (`AGENT_PLAN.md` section 3): the library
/// returns the typed [`methyltfr::Error`] and only the binary flattens it for
/// display. A non-zero exit is what a drop-in replacement must do on the inputs
/// upstream aborts on.
fn main() -> Result<()> {
    main_with(Cli::parse()).map_err(|e| anyhow::anyhow!("{e}"))
}
