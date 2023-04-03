use clap::Parser;
use motif_finder::{Error, MotifFinder};
use tracing::Level;
use tracing_subscriber::FmtSubscriber;

fn main() -> Result<(), Error> {


    

    let motif_finder = MotifFinder::parse();
    let subscriber = FmtSubscriber::builder()
    .pretty()
    .with_max_level(Level::TRACE)
    .with_target(false)
    .with_line_number(false)
    .compact()
    .with_file(false)
    .without_time()
    .finish();
    tracing::subscriber::set_global_default(subscriber).expect("setting default subscriber failed");
    motif_finder.exec()
}
