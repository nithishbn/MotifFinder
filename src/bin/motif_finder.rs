use clap::Parser;

use motif_finder::{Error, MotifFinder};

fn main() -> Result<(), Error> {
    let motif_finder = MotifFinder::parse();

    motif_finder.exec()
}
