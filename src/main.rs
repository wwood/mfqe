use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::fs::File;
use std::fs::OpenOptions;
use std::collections::{HashMap,HashSet};
use std::env;
use std::path::Path;

extern crate seq_io;
use seq_io::fastq::Record;
use seq_io::fasta::Record as OtherRecord;

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
    let app = Command::new("mfqe")
        .version(crate_version!())
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Extract multiple sets of fastq reads by name")
        .override_usage("\nUsage for FASTQ:\n  \
                  zcat my.fastq.gz |mfqe --sequence-name-lists <LIST1> .. --output-fastq-files <OUTPUT1> ..\n\
               \n\
\n\
               Extract one or more sets of reads from a FASTQ (or \
               FASTA) file by specifying their read names.\n\n\
\n\
               Read name files are uncompressed text files with read names \
               (without comments).\n\
\n\
               Output is gzip-compressed unless --output-uncompressed is specified, input is uncompressed.\n\
               \nOther FASTQ options:
               \n--input-fastq <PATH>: Use this file as input FASTQ [default: Use STDIN]\
               \n\n\
               An analogous set of options is implemented for FASTA:\n\n\
               --output-fasta-files <OUTPUT1> ..\n\
               --input-fasta <PATH>\n\n")

        // Unfortunately clap cannot properly take multiple .long() arguments, so
        // we have to define multiple args for backwards compatibility.
        .arg(Arg::new("sequence-name-lists")
             .long("sequence-name-lists")
             .short('l')
             .help("List of files each containing sequence IDs")
             .required_unless_present_any(&["fastq-read-name-lists","fasta-read-name-lists"])
             .conflicts_with_all(&["fastq-read-name-lists","fasta-read-name-lists"])
             .takes_value(true)
             .multiple_occurrences(true)
             .multiple_values(true))
        .arg(Arg::new("fastq-read-name-lists")
             .long("fastq-read-name-lists")
             .help("List of files each containing sequence IDs [alias for --sequence-name-lists]")
             .required_unless_present_any(&["sequence-name-lists","fasta-read-name-lists"])
             .conflicts_with_all(&["sequence-name-lists","fasta-read-name-lists"])
             .takes_value(true)
             .multiple_occurrences(true)
             .multiple_values(true))
        .arg(Arg::new("fasta-read-name-lists")
             .long("fasta-read-name-lists")
             .help("List of files each containing sequence IDs [alias for --sequence-name-lists]")
             .required_unless_present_any(&["fastq-read-name-lists","sequence-name-lists"])
             .conflicts_with_all(&["fastq-read-name-lists","sequence-name-lists"])
             .takes_value(true)
             .multiple_occurrences(true)
             .multiple_values(true))

        .arg(Arg::new("output-fastq-files")
             .long("output-fastq-files")
             .help("List of files to write FASTQ to")
             .required_unless_present("output-fasta-files")
             .takes_value(true)
             .multiple_occurrences(true)
             .multiple_values(true))
        .arg(Arg::new("input-fastq")
             .long("input-fastq")
             .help("File containing uncompressed input FASTQ sequences [default: Use STDIN]")
             .takes_value(true))

        .arg(Arg::new("output-fasta-files")
             .long("output-fasta-files")
             .help("List of files to write FASTA to")
             .required_unless_present("output-fastq-files")
             .takes_value(true)
             .multiple_occurrences(true)
             .multiple_values(true))
        .arg(Arg::new("input-fasta")
             .long("input-fasta")
             .help("File containing uncompressed input FASTA sequences [default: Use STDIN]")
             .takes_value(true))
             
        .arg(Arg::new("output-uncompressed")
             .long("output-uncompressed")
             .help("Output sequences uncompressed [default: gzip compress outputs]")
             .short('u'))
        .arg(Arg::new("append")
             .long("append")
             .help("Append to output files [default: Overwrite]")
             .short('a'));

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

    let output_files: Vec<&str>;
    let input: Option<BufReader<File>>;

    let doing_fastq = matches.is_present("output-fastq-files");
    let read_lists: Vec<&str> = if matches.is_present("fasta-read-name-lists") {
        matches.values_of("fasta-read-name-lists").unwrap().collect()
    } else if matches.is_present("fastq-read-name-lists") {
        matches.values_of("fastq-read-name-lists").unwrap().collect()
    } else {
        matches.values_of("sequence-name-lists").unwrap().collect()
    };
    debug!("Found readname lists {:#?}", read_lists);

    let appending = matches.is_present("append");

    if doing_fastq {
        // Doing fastq
        output_files = matches.values_of("output-fastq-files").unwrap().collect();
        input = match matches.value_of("input-fastq") {
            Some(path) => Some(BufReader::new(
                File::open(path)
                    .expect("Failed to open fastq file for reading"))),
            None => None
        };

    } else {
        output_files = matches.values_of("output-fasta-files").unwrap().collect();
        input = match matches.value_of("input-fasta") {
            Some(path) => Some(BufReader::new(
                File::open(path)
                    .expect("Failed to open fasta file for reading"))),
            None => None
        };

    }

    // check the number of fastq read name files is the same as the number of
    // output files.
    if output_files.len() != read_lists.len() {
        panic!("The number of read name files was {}, output files there was \
                {}. These numbers must be equal.",
               output_files.len(), read_lists.len());
    }

    let name_index = generate_name_index(read_lists);

    // Open output file as gzipped output
    let output_compressed = !matches.is_present("output-uncompressed");
    let uncompressed_outputs: Option<Vec<File>> = if output_compressed {
        None 
    } else {
        Some(output_files.iter().map( |o| {
            open_a_file(o, appending)
        }).collect())
    };
    let compressed_outputs: Option<Vec<GzEncoder<BufWriter<File>>>> = if output_compressed {
        Some(output_files.iter().map( |o| {
            let w1 = open_a_file(o, appending);
            GzEncoder::new(BufWriter::new(w1), Compression::default())
            }).collect())
    } else {
        None
    };

    match doing_fastq {
        true => match output_compressed {
            true => fastq_pipeline(input, name_index, compressed_outputs.unwrap()),
            false => fastq_pipeline(input, name_index, uncompressed_outputs.unwrap()),
        },
        false => match output_compressed {
            true => fasta_pipeline(input, name_index, compressed_outputs.unwrap()),
            false => fasta_pipeline(input, name_index, uncompressed_outputs.unwrap()),
        },
    };
}

fn open_a_file(path: &str, append: bool) -> File {
    match append {
        true => {
            // must create the file, if it doesn't already exist
            if !Path::new(path).exists() {
                File::create(path).expect(&format!("Failed to create output file {}, before opening in append mode", path));
            }
            OpenOptions::new().append(true).open(path)
                .expect(&format!("Failed to open output file {} for appending", path))
        },
        false => File::create(path)
            .expect(&format!("Failed to open output file {} for writing", path)),
    }
}

struct NameIndex {
    name_to_index: HashMap<String, HashSet<usize>>,
    index_to_expected_count: Vec<usize>,
}

fn generate_name_index(read_lists: Vec<&str>) -> NameIndex {
    // Read in each read name into has hashmap
    let mut name_to_index: HashMap<String, HashSet<usize>> = HashMap::new();
    let mut index_to_expected_count: Vec<usize> = vec![];
    for (i, read_name_file) in read_lists.iter().enumerate() {
        let mut lines_in_file: u64 = 0;
        let reader1 = File::open(read_name_file)
            .expect(&format!("Failed to open read name file {}", read_name_file));
        let reader = BufReader::new(reader1);
        for line in reader.lines() {
            let name = line.unwrap();
            // Ignore blank lines
            if name != "" {
                let insert = match name_to_index.get_mut(&name) {
                    Some(prevs) => {
                        if !prevs.insert(i) {
                            panic!(
                                "It appears that read '{}' was specified twice in input file {}",
                                name, read_lists[i]);
                        }
                        false
                    },
                    None => true
                };
                if insert { // Do this to get around the borrow checker
                    let mut set = HashSet::with_capacity(1);
                    set.insert(i);
                    name_to_index.insert(name.clone(), set);
                }
                lines_in_file += 1;
            }
        }
        index_to_expected_count.push(lines_in_file as usize);
        info!("Read in {} read names from {}", lines_in_file, read_name_file);
    }

    NameIndex {
        name_to_index: name_to_index,
        index_to_expected_count: index_to_expected_count
    }
}
fn fastq_pipeline<W: Write>(
    fastq_input: Option<BufReader<File>>,
    name_index: NameIndex,
    outputs: Vec<W>) {

    match fastq_input {
        Some(r) => read_fastq(
            seq_io::fastq::Reader::new(r),
            name_index.index_to_expected_count,
            name_index.name_to_index,
            outputs),
        None => read_fastq(
            seq_io::fastq::Reader::new(std::io::stdin()),
            name_index.index_to_expected_count,
            name_index.name_to_index,
            outputs)
    };
}

fn read_fastq<R, W>(
    mut reader: seq_io::fastq::Reader<R>,
    index_to_expected_count: Vec<usize>,
    name_to_index: HashMap<String, HashSet<usize>>,
    mut fastq_outputs: Vec<W>)
where R: Read, W: Write {
    info!("Iterating input FASTQ file");
    let mut total_input_reads: usize = 0;
    let mut index_to_observed_count: Vec<usize> = vec![0; index_to_expected_count.len()];

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

fn fasta_pipeline<W: Write>(
    input: Option<BufReader<File>>,
    name_index: NameIndex,
    outputs: Vec<W>) {


    match input {
        Some(r) => read_fasta(
            seq_io::fasta::Reader::new(r),
            name_index.index_to_expected_count,
            name_index.name_to_index,
            outputs),
        None => read_fasta(
            seq_io::fasta::Reader::new(std::io::stdin()),
            name_index.index_to_expected_count,
            name_index.name_to_index,
            outputs)
    };

}


fn read_fasta<R, W>( // TODO: This is duplicated code, but too lazy to fix right now.
    mut reader: seq_io::fasta::Reader<R>,
    index_to_expected_count: Vec<usize>,
    name_to_index: HashMap<String, HashSet<usize>>,
    mut fastq_outputs: Vec<W>)
where R: Read, W: Write {
    info!("Iterating input FASTQ file");
    let mut total_input_reads: usize = 0;
    let mut index_to_observed_count: Vec<usize> = vec![0; index_to_expected_count.len()];

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
