use std::fs::{File, self};
use std::io::{self, BufRead, Write};
use std::path::Path;

use chrono::Utc;
use clap::Parser;
use Motif_Finder::{iterate_gibbs_sampler, median_string::median_string};

/// Motif Finder
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// file path of the genome
    #[arg(short, long = "path")]
    path_to_file: String,

    /// how many entries to read
    #[arg(short = 'e', long = "entry")]
    num_entries: usize,

    /// motif length
    #[arg(short)]
    k: usize,

    /// number of runs
    #[arg(short = 'r', long = "runs", default_value_t = 10)]
    num_runs: usize,

    /// number of iterations per run
    #[arg(short = 'i', long = "iters", default_value_t = 10)]
    num_iterations: usize,

    /// save motifs to file
    #[arg(short)]
    save: bool,
}

fn main() {
    let args = Args::parse();
    let dt = Utc::now();
    let timestamp: i64 = dt.timestamp();
    println!("start at {}", dt);
    // identify_ori("./p_tricornutum.fa",1);
    let mut genome = vec![];
    let mut chromosome = String::from("");
    let mut current = 1;
    println!("Loading data...");
    if let Ok(lines) = read_lines(args.path_to_file) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(line) = line {
                if line.starts_with(">") {
                    if chromosome.chars().count() > 0 {
                        // dbg!("hi");
                        // dbg!(line);
                        if chromosome.chars().count() == 0{
                            continue;
                        }
                        genome.push(chromosome);
                        chromosome = "".to_string();
                        current += 1;
                        continue;
                    }else if chromosome.chars().count() == 0{
                        continue;
                    }
                }
                if current == args.num_entries + 1{
                    if chromosome.chars().count() == 0{
                        continue;
                    }
                    genome.push(chromosome);
                    chromosome = "".to_string();
                    break;
                }
                
                chromosome += line.to_uppercase().trim();
            }
        }
    }

    
    if chromosome.chars().count() > 0 {
        genome.push(chromosome);
    }
    println!("done loading data: {} chromosomes", genome.len());
    
    // let median_string = median_string(args.k, &genome, timestamp);
    // println!("median string: {}", median_string);
    let motifs = iterate_gibbs_sampler(
        &genome,
        args.k,
        genome.len(),
        args.num_iterations,
        args.num_runs,
    );
    if args.save {
        let mut file = fs::File::create(format!("MotifFinder-output-{timestamp}-{}.txt", args.k))
            .expect("Unable to create file");
        let version = env!("CARGO_PKG_VERSION");
        writeln!(file, "MotifFinder {}",version).expect("Unable to write data");
        writeln!(file, "k: {}", args.k).expect("Unable to write data");
        writeln!(file, "number of entries: {}", args.num_entries).expect("Unable to write data");
        writeln!(file, "runs: {}", args.num_runs).expect("Unable to write data");
        writeln!(file, "iterations: {}", args.num_iterations).expect("Unable to write data");
        writeln!(file, "Start time: {}", dt).expect("Unable to write data");
        
        
        let dt = Utc::now();
        writeln!(file, "End time: {}", dt).expect("Unable to write data");
        writeln!(file,"_______________________________________________________________________________________").expect("Unable to write data");

        for (i,motif) in motifs.iter().enumerate(){
            let motif = motif.trim();
            write!(file, ">motif {}\n", i+1).expect("Unable to write data");
            write!(file, "{}\n", motif).expect("Unable to write data");
        }
    }
    for motif in motifs {
        println!("{}", motif);
    }
    let dt = Utc::now();
    println!("done at {}", dt.timestamp() - timestamp);
    
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}