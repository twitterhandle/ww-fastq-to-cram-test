version 1.0
# Convert a group of paired fastq.gz files into an unmapped cram
#  Uses the convention: READ_GROUP_NAME=~{sampleName}_~{flowcellName}

#### STRUCT DEFINITIONS

struct flowcellFastqs {
  String flowcellName # name of flowcell 
  Array[File] fastqR1Locations # array of input R1 fastq file locations
  Array[File] fastqR2Locations # array of input R2 fastq file locations
}

struct inputData {
  String datasetID # unique ID to identify a particular dataset, even if sampleName is not unique
  String sampleName # sample name to insert into the read group header
  String libraryName # library name to place into the LB attribute in the read group header
  String sequencingCenter # location where the sample was sequenced
  Array[flowcellFastqs] filepaths # array of flowcell fastq file locations
}

#### WORKFLOW DEFINITION

workflow PairedFastqsToUnmappedCram {
  input {
    Array[inputData] batchInfo
  }

  String GATKDocker = "ghcr.io/getwilds/gatk:4.3.0.0"
  String samtoolsDocker = "ghcr.io/getwilds/samtools:1.11"

  scatter (sample in batchInfo) { # for every sample in your batch,
    String base_file_name = sample.sampleName + "_" + sample.datasetID

    scatter (flowcell in sample.filepaths) { # and for every flowcell that sample's library was run on,  
      call FastqtoUnmappedBam { # take all the fastqs for that sample in that flowcell and make an unmapped bam
        input:
          R1fastq = flowcell.fastqR1Locations,
          R2fastq = flowcell.fastqR2Locations,
          base_file_name = base_file_name + "_" + flowcell.flowcellName,
          sampleName = sample.sampleName,
          sequencingCenter = sample.sequencingCenter,
          libraryName = sample.libraryName,
          flowcellName = flowcell.flowcellName,
          docker = GATKDocker
      }
    } # End flowcell scatter

    call mergeBamstoCram { # then for each flowcell, merge all unmapped bams into one unmapped cram for the library
      input:
        bamsToMerge = FastqtoUnmappedBam.unmappedbam,
        base_file_name = base_file_name,
        docker = samtoolsDocker,
        threads = 6
    }

    call ValidateCram { # then validate to make sure it checks out
      input: 
        unmappedCram=mergeBamstoCram.cram,
        base_file_name = base_file_name,
        docker = GATKDocker
    }
  } # End sample scatter

  # Outputs that will be retained when execution is complete
  output {
    Array[File] unmappedCrams = mergeBamstoCram.cram
    Array[File] unmappedCramIndexes = mergeBamstoCram.crai
    Array[File] validation = ValidateCram.validation
  }

  parameter_meta (
    batchInfo: "array of inputData structs describing the relevant metadata for each sample"

    unmappedCrams: "array of unmapped cram files for each sample"
    unmappedCramIndexes: "array of index files for each unmapped cram file"
    validation: "text file containing all relevant validation statistics for the cram in question"
  )
} # End workflow

#### TASK DEFINITIONS

# Converts fastq files to unaligned bams
task FastqtoUnmappedBam {
  input {
    Array[File] R1fastq
    Array[File] R2fastq
    String base_file_name
    String sampleName
    String flowcellName
    String libraryName
    String sequencingCenter
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
    FastqToSam \
      --FASTQ ~{sep=" " R1fastq} \
      --FASTQ2 ~{sep=" " R2fastq} \
      --OUTPUT ~{base_file_name}.unmapped.bam \
      --READ_GROUP_NAME ~{sampleName}_~{flowcellName} \
      --SAMPLE_NAME ~{sampleName} \
      --LIBRARY_NAME ~{libraryName} \
      --PLATFORM illumina \
      --SEQUENCING_CENTER ~{sequencingCenter}
  >>>

  output {
    File unmappedbam = "~{base_file_name}.unmapped.bam"
  }

  runtime {
    cpu: 4
    memory: "8G"
    docker: docker
  }

  parameter_meta {
    R1fastq: "array of R1 fastq files for the library in question"
    R2fastq: "array of R2 fastq files for the library in question"
    base_file_name: "base file name to use in the read group and bam file names"
    sampleName: "name of the sample in question"
    flowcellName: "name of the flowcell on which the sample is being sequenced"
    libraryName: "name of the library in question"
    sequencingCenter: "location where the sample was sequenced"
    docker: "location of Docker image to use for this task"

    unmappedbam: "final unmapped bam file containing the reads from the fastqs in question"
  }
}

# Validates cram files for formatting issues. 
task ValidateCram {
  input {
    File unmappedCram
    String base_file_name
    String docker
  }

  command <<<
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms2g" \
      ValidateSamFile \
        --INPUT ~{unmappedCram} \
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
    unmappedCram: "unmapped cram file to validate"
    base_file_name: "base file name to use when saving the validation text file"
    docker: "location of Docker image to use for this task"

    validation: "text file containing all relevant validation statistics for the cram in question"
  }
}

# Merges multiple bam files into a single cram file
task mergeBamstoCram {
  input {
    Array[File] bamsToMerge
    String base_file_name
    String docker
    Int threads
  }

  command <<<
    set -eo pipefail
    samtools merge -@ ~{threads-1} \
      --write-index --output-fmt CRAM  \
      ~{base_file_name}.merged.cram ~{sep=" " bamsToMerge} 
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
    bamsToMerge: "array of bam files to merge into a single cram file"
    base_file_name: "base file name to use when saving the cram file"
    docker: "location of Docker image to use for this task"
    threads: "number of threads to use during the merging process"

    cram: "final cram file containing all reads"
    crai: "index file of final cram"
  }
}