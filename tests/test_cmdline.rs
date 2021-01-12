extern crate assert_cli;
extern crate tempfile;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    extern crate tempfile;
    use std::io::Read;

    #[test]
    fn test_fastq_by_file(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fastq-read-name-lists",
                "tests/data/input1",
                "--output-fastq-files",
                t,
                "--input-fastq",
                "tests/data/1.fq"]).succeeds().unwrap();
        Assert::command(&["zcat",t])
            .stdout().is("@random_sequence_length_5_1 1\n\
                          TAGGG\n\
                          +\n\
                          AAAAA\n").unwrap();
    }

    #[test]
    fn test_fastq_by_stdin(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        let mut contents = String::new();
        std::fs::File::open("tests/data/1.fq").unwrap().read_to_string(&mut contents).unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fastq-read-name-lists",
                "tests/data/input2",
                "--output-fastq-files",
                t])
            .stdin(contents)
            .succeeds().unwrap();
        Assert::command(&["zcat",t])
            .stdout().is("@random_sequence_length_5_1 1\n\
                          TAGGG\n\
                          +\n\
                          AAAAA\n\
                          @random_sequence_length_5_2 2\n\
                          TTTCA\n\
                          +\n\
                          ATGCA\n").unwrap();
    }

    #[test]
    fn test_fasta_by_file(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fasta-read-name-lists",
                "tests/data/input1",
                "--output-fasta-files",
                t,
                "--input-fasta",
                "tests/data/a.fasta"]).succeeds().unwrap();
        Assert::command(&["zcat",t])
            .stdout().is(">random_sequence_length_5_1\n\
                          GGTGT\n").unwrap();
    }

    #[test]
    fn test_fasta_by_stdin(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        let mut contents = String::new();
        std::fs::File::open("tests/data/a.fasta").unwrap().read_to_string(&mut contents).unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fasta-read-name-lists",
                "tests/data/input1",
                "--output-fasta-files",
                t])
            .stdin(contents)
            .succeeds().unwrap();
        Assert::command(&["zcat",t])
            .stdout().is(">random_sequence_length_5_1\n\
                          GGTGT\n").unwrap();
    }

    #[test]
    fn test_define_sequence_names(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        let mut contents = String::new();
        std::fs::File::open("tests/data/a.fasta").unwrap().read_to_string(&mut contents).unwrap();
        Assert::main_binary()
            .with_args(&[
                "-l",
                "tests/data/input1",
                "--output-fasta-files",
                t])
            .stdin(contents)
            .succeeds().unwrap();
        Assert::command(&["zcat",t])
            .stdout().is(">random_sequence_length_5_1\n\
                          GGTGT\n").unwrap();
    }

    #[test]
    fn test_fastq_by_file_empty_lines_in_read_names(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fastq-read-name-lists",
                "tests/data/input1_with_empty",
                "--output-fastq-files",
                t,
                "--input-fastq",
                "tests/data/1.fq"]).succeeds().unwrap();
        Assert::command(&["zcat",t])
            .stdout().is("@random_sequence_length_5_1 1\n\
                          TAGGG\n\
                          +\n\
                          AAAAA\n").unwrap();
    }

}
