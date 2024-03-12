version 1.0
# Convert a group of paired fastq.gz files into an unmapped cram
#  Uses the convention: READ_GROUP_NAME=~{sample_name}_~{flowcell_name}

#### STRUCT DEFINITIONS

struct FlowcellFastqs {
  String flowcell_name # name of flowcell 
  Array[File] fastq_r1_locations # array of input R1 fastq file locations
  Array[File] fastq_r2_locations # array of input R2 fastq file locations
}

struct InputData {
  String dataset_id # unique ID to identify a particular dataset, even if sample_name is not unique
  String sample_name # sample name to insert into the read group header
  String library_name # library name to place into the LB attribute in the read group header
  String sequencing_center # location where the sample was sequenced
  Array[FlowcellFastqs] filepaths # array of flowcell fastq file locations
}

#### WORKFLOW DEFINITION

workflow PairedFastqsToUnmappedCram {
  input {
    Array[InputData] batch_info
  }

  String gatk_docker = "ghcr.io/getwilds/gatk:4.3.0.0"
  String samtools_docker = "ghcr.io/getwilds/samtools:1.11"

  scatter (sample in batch_info) { # for every sample in your batch,
    String base_file_name = sample.sample_name + "_" + sample.dataset_id

    scatter (flowcell in sample.filepaths) { # and for every flowcell that sample's library was run on,  
      call FastqToUnmappedBam { # take all the fastqs for that sample in that flowcell and make an unmapped bam
        input:
          r1_fastq = flowcell.fastq_r1_locations,
          r2_fastq = flowcell.fastq_r2_locations,
          base_file_name = base_file_name + "_" + flowcell.flowcell_name,
          sample_name = sample.sample_name,
          sequencing_center = sample.sequencing_center,
          library_name = sample.library_name,
          flowcell_name = flowcell.flowcell_name,
          docker = gatk_docker
      }
    } # End flowcell scatter

    call MergeBamsToCram { # then for each flowcell, merge all unmapped bams into one unmapped cram for the library
      input:
        bams_to_merge = FastqToUnmappedBam.unmapped_bam,
        base_file_name = base_file_name,
        docker = samtools_docker,
        threads = 6
    }

    call ValidateCram { # then validate to make sure it checks out
      input: 
        unmapped_cram=MergeBamsToCram.cram,
        base_file_name = base_file_name,
        docker = gatk_docker
    }
  } # End sample scatter

  # Outputs that will be retained when execution is complete
  output {
    Array[File] unmapped_crams = MergeBamsToCram.cram
    Array[File] unmapped_cram_indexes = MergeBamsToCram.crai
    Array[File] validation = ValidateCram.validation
  }

  parameter_meta {
    batch_info: "array of InputData structs describing the relevant metadata for each sample"

    unmapped_crams: "array of unmapped cram files for each sample"
    unmapped_cram_indexes: "array of index files for each unmapped cram file"
    validation: "text file containing all relevant validation statistics for the cram in question"
  }
} # End workflow

#### TASK DEFINITIONS

# Converts fastq files to unaligned bams
task FastqToUnmappedBam {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    String base_file_name
    String sample_name
    String flowcell_name
    String library_name
    String sequencing_center
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
    FastqToSam \
      --FASTQ ~{sep=" " r1_fastq} \
      --FASTQ2 ~{sep=" " r2_fastq} \
      --OUTPUT ~{base_file_name}.unmapped.bam \
      --READ_GROUP_NAME ~{sample_name}_~{flowcell_name} \
      --SAMPLE_NAME ~{sample_name} \
      --LIBRARY_NAME ~{library_name} \
      --PLATFORM illumina \
      --SEQUENCING_CENTER ~{sequencing_center}
  >>>

  output {
    File unmapped_bam = "~{base_file_name}.unmapped.bam"
  }

  runtime {
    cpu: 4
    memory: "8G"
    docker: docker
  }

  parameter_meta {
    r1_fastq: "array of R1 fastq files for the library in question"
    r2_fastq: "array of R2 fastq files for the library in question"
    base_file_name: "base file name to use in the read group and bam file names"
    sample_name: "name of the sample in question"
    flowcell_name: "name of the flowcell on which the sample is being sequenced"
    library_name: "name of the library in question"
    sequencing_center: "location where the sample was sequenced"
    docker: "location of Docker image to use for this task"

    unmapped_bam: "final unmapped bam file containing the reads from the fastqs in question"
  }
}

# Validates cram files for formatting issues. 
task ValidateCram {
  input {
    File unmapped_cram
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms2g" \
      ValidateSamFile \
        --INPUT ~{unmapped_cram} \
        --MODE SUMMARY \
        --IGNORE_WARNINGS false > ~{base_file_name}.validation.txt
  >>>

  output {
    File validation = "~{base_file_name}.validation.txt"
  }

  runtime {
    cpu: 2
    memory: "4 GB"
    docker: docker
  }

  parameter_meta {
    unmapped_cram: "unmapped cram file to validate"
    base_file_name: "base file name to use when saving the validation text file"
    docker: "location of Docker image to use for this task"

    validation: "text file containing all relevant validation statistics for the cram in question"
  }
}

# Merges multiple bam files into a single cram file
task MergeBamsToCram {
  input {
    Array[File] bams_to_merge
    String base_file_name
    String docker
    Int threads
  }

  command <<<
    set -eo pipefail
    samtools merge -@ ~{threads-1} \
      --write-index --output-fmt CRAM  \
      ~{base_file_name}.merged.cram ~{sep=" " bams_to_merge} 
  >>>

  output {
    File cram = "~{base_file_name}.merged.cram"
    File crai = "~{base_file_name}.merged.cram.crai"
  }

  runtime {
    docker: docker
    cpu: threads
  }

  parameter_meta {
    bams_to_merge: "array of bam files to merge into a single cram file"
    base_file_name: "base file name to use when saving the cram file"
    docker: "location of Docker image to use for this task"
    threads: "number of threads to use during the merging process"

    cram: "final cram file containing all reads"
    crai: "index file of final cram"
  }
}