use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::fs::File;
use std::collections::{HashMap,HashSet};
use std::env;

extern crate seq_io;
use seq_io::fastq::Record;

extern crate clap;
use clap::*;

extern crate flate2;
use flate2::Compression;
use flate2::write::GzEncoder;

#[macro_use]
extern crate log;
extern crate env_logger;
use log::LevelFilter;
use env_logger::Builder;


fn main() {
    let app = App::new("mfqe")
        .version(crate_version!())
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Extract multiple sets of fastq reads by name")
        .help("Usage: zcat my.fastq.gz |mfqe --fastq-read-name-lists <LIST1> .. --output-fastq-files <OUTPUT1> ..\n\
               \n\
               Purpose is to extract one or more sets of reads from a FASTQ file by specifying their read name.\n\n\
               Read name files are uncompressed text files with read names (without comments).\n\
               Output is gzip-compressed, input may or may not be.\n")

        .arg(Arg::with_name("fastq-read-name-lists")
             .long("fastq-read-name-lists")
             .required(true)
             .takes_value(true)
             .multiple(true))
        .arg(Arg::with_name("output-fastq-files")
             .long("output-fastq-files")
             .required(true)
             .takes_value(true)
             .multiple(true));

    let matches = app.clone().get_matches();

    let log_level = LevelFilter::Info;
    let mut builder = Builder::new();
    builder.filter_level(log_level);
    if env::var("RUST_LOG").is_ok() {
        builder.parse_filters(&env::var("RUST_LOG").unwrap());
    }
    if builder.try_init().is_err() {
        panic!("Failed to set log level - has it been specified multiple times?")
    }

    // check the number of fastq read name files is the same as the number of
    // output files.
    let fastq_read_lists: Vec<&str> = matches.values_of("fastq-read-name-lists").unwrap().collect();
    let fastq_output_files: Vec<&str> = matches.values_of("output-fastq-files").unwrap().collect();

    if fastq_output_files.len() != fastq_read_lists.len() {
        panic!("The number of read name files was {}, output files there was \
                {}. These numbers must be equal.",
               fastq_output_files.len(), fastq_read_lists.len());
    }

    // Read in each read name into has hashmap
    let mut name_to_index: HashMap<String, HashSet<usize>> = HashMap::new();
    let mut index_to_expected_count: Vec<usize> = vec![];
    for (i, read_name_file) in fastq_read_lists.iter().enumerate() {
        let mut lines_in_file: u64 = 0;
        let reader1 = File::open(read_name_file)
            .expect(&format!("Failed to open read name file {}", read_name_file));
        let reader = BufReader::new(reader1);
        for line in reader.lines() {
            let name = line.unwrap();
            match name_to_index.get_mut(&name) {
                Some(prevs) => {
                    if !prevs.insert(i) {
                        panic!(
                            "It appears that read '{}' was specified twice in input file {}",
                            name, fastq_read_lists[i]);
                    }
                },
                None => {
                    let mut set = HashSet::with_capacity(1);
                    set.insert(i);
                    name_to_index.insert(name.clone(), set);
                }
            }
            lines_in_file += 1;
        }
        index_to_expected_count.push(lines_in_file as usize);
        info!("Read in {} read names from {}", lines_in_file, read_name_file);
    }

    // Open output file as gzipped output
    info!("Opening output FASTQ files ..");
    let mut fastq_outputs: Vec<GzEncoder<BufWriter<File>>> = fastq_output_files.iter().map( |o| {
        let w1 = File::create(o)
            .expect(&format!("Failed to open output file {} for writing", o));
        GzEncoder::new(BufWriter::new(w1), Compression::default())
    }).collect();

    info!("Iterating input FASTQ file");
    let mut total_input_reads: usize = 0;
    let mut index_to_observed_count: Vec<usize> = vec![0; index_to_expected_count.len()];
    let mut reader = seq_io::fastq::Reader::new(std::io::stdin());
    while let Some(record) = reader.next() {
        let r2 = record.unwrap();
        match name_to_index.get(r2.id()
                                .expect("UTF8 error when decoding FASTQ header")) {
            Some(indices) => {
                for i in indices {
                    index_to_observed_count[*i] += 1;
                    r2.write(&mut fastq_outputs[*i]).expect("Failed to write a FASTQ record");
                }
            },
            None => {}
        };
        total_input_reads += 1;
    }

    let total_assigned_reads: usize = index_to_observed_count.iter().sum();
    info!("Extracted {} reads from {} total", total_assigned_reads, total_input_reads);
    if index_to_expected_count != index_to_observed_count {
        panic!("Mismatching numbers of read names were observed. Expected:\n{:?}\nbut found\n{:?}",
               index_to_expected_count, index_to_observed_count);
    }
}
