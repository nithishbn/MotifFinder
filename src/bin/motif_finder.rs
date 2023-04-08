use clap::Parser;
use motif_finder::{Error, MotifFinder};
use tracing::Level;
use tracing_subscriber::FmtSubscriber;

fn main() -> Result<(), Error> {
    let motif_finder = MotifFinder::parse();
    let log_level = motif_finder
        .verbose
        .log_level()
        .ok_or(Error::GenericError)?;
    let converted_level = match log_level {
        clap_verbosity_flag::Level::Error => Level::ERROR,
        clap_verbosity_flag::Level::Warn => Level::WARN,
        clap_verbosity_flag::Level::Info => Level::INFO,
        clap_verbosity_flag::Level::Debug => Level::DEBUG,
        clap_verbosity_flag::Level::Trace => Level::TRACE,
    };
    let subscriber = FmtSubscriber::builder()
        .pretty()
        .with_max_level(converted_level)
        .with_target(false)
        .with_line_number(false)
        .compact()
        .with_file(false)
        .without_time()
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("setting default subscriber failed");
    motif_finder.exec()
}
