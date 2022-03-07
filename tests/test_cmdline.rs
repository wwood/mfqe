extern crate assert_cli;
extern crate tempfile;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    extern crate tempfile;
    use std::io::Read;
    use std::io::Write;

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

    #[test]
    fn test_appending_gzip(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        let mut contents = String::new();
        std::fs::File::open("tests/data/a.fasta").unwrap().read_to_string(&mut contents).unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fasta-read-name-lists",
                "tests/data/input1",
                "--output-fasta-files",
                t,
                "--append",
                "--input-fasta",
                "tests/data/a.fasta"]).succeeds().unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fasta-read-name-lists",
                "tests/data/input1",
                "--output-fasta-files",
                t,
                "--append",
                "--input-fasta",
                "tests/data/a.fasta"]).succeeds().unwrap();
        Assert::command(&["zcat",t])
            .stdout().is(">random_sequence_length_5_1\n\
            GGTGT\n\
            >random_sequence_length_5_1\n\
            GGTGT\n").unwrap();
    }


    #[test]
    fn test_appending_no_gzip(){
        let mut tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        write!(tf, "abc\n").unwrap();
        tf.flush().unwrap();
        let t = tf.path().to_str().unwrap();
        let mut contents = String::new();
        std::fs::File::open("tests/data/a.fasta").unwrap().read_to_string(&mut contents).unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fasta-read-name-lists",
                "tests/data/input1",
                "--output-fasta-files",
                t,
                "--append",
                "--output-uncompressed",
                "--input-fasta",
                "tests/data/a.fasta"]).succeeds().unwrap();
        Assert::command(&["cat",t])
            .stdout().is("abc\n>random_sequence_length_5_1\n\
                          GGTGT\n").unwrap();
    }

    #[test]
    fn test_appending_not_already_existing(){
        let t = "/tmp/testme123_mfqe";
        // delete t if it already exists
        if std::path::Path::new(t).exists() {
            std::fs::remove_file(t).unwrap();
        }
        let mut contents = String::new();
        std::fs::File::open("tests/data/a.fasta").unwrap().read_to_string(&mut contents).unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fasta-read-name-lists",
                "tests/data/input1",
                "--output-fasta-files",
                t,
                "--append",
                "--output-uncompressed",
                "--input-fasta",
                "tests/data/a.fasta"]).succeeds().unwrap();
        Assert::command(&["cat",t])
            .stdout().is(">random_sequence_length_5_1\n\
                          GGTGT\n").unwrap();
    }


    #[test]
    fn test_appending_no_gzip_two_files(){
        let mut tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        write!(tf, "abc\n").unwrap();
        tf.flush().unwrap();

        let mut tf2: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        write!(tf2, "defabc\n").unwrap();
        tf2.flush().unwrap();
        
        let t = tf.path().to_str().unwrap();
        let t2 = tf2.path().to_str().unwrap();

        let mut contents = String::new();
        std::fs::File::open("tests/data/a.fasta").unwrap().read_to_string(&mut contents).unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fasta-read-name-lists",
                "tests/data/input1",
                "tests/data/input2",
                "--output-fasta-files",
                t,
                t2,
                "--append",
                "--output-uncompressed",
                "--input-fasta",
                "tests/data/a.fasta"]).succeeds().unwrap();
        Assert::command(&["cat",t])
            .stdout().is("abc\n>random_sequence_length_5_1\n\
                            GGTGT\n").succeeds().unwrap();
        Assert::command(&["cat",t2])
            .stdout().is("defabc\n>random_sequence_length_5_1\n\
            GGTGT\n\
            >random_sequence_length_5_2\n\
            TTATG\n").succeeds().unwrap();
    }


    #[test]
    fn test_fasta_sequence_prefix(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        let mut contents = String::new();
        std::fs::File::open("tests/data/a.fasta").unwrap().read_to_string(&mut contents).unwrap();
        Assert::main_binary()
            .with_args(&[
                "--fasta-read-name-lists",
                "tests/data/input1",
                "--sequence-prefix",
                "i am a prefix",
                "--output-fasta-files",
                t])
            .stdin(contents)
            .succeeds().unwrap();
        Assert::command(&["zcat",t])
            .stdout().is(">i am a prefixrandom_sequence_length_5_1\n\
                          GGTGT\n").unwrap();
    }
}
