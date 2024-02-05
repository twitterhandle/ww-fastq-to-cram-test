version 1.0
# Convert a group of paired fastq.gz files into an unmapped cram
#  Uses the convention: READ_GROUP_NAME=~{sampleName}_~{flowcellName}

struct flowcellFastqs {
  String flowcellName
  Array[File] fastqR1Locations
  Array[File] fastqR2Locations
}

struct inputData {
  String datasetID # a unique ID for this data set that would identify a particular dataset even if smapleName was not unique.
  String sampleName # the sample name  - this is often what is used as the beginning of your filenames
  String libraryName
  String sequencingCenter
  Array[flowcellFastqs] filepaths
}

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
    call mergeBamstoCram { # then for each of the flowcells that library was run on, merge all the unmapped bams into one unmapped bam for the library
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
# End workflow
}

#### TASK DEFINITIONS

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
  command {
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
  }
  runtime {
    cpu: 4
    memory: "8G"
    docker: docker
  }
  output {
    File unmappedbam = "~{base_file_name}.unmapped.bam"
  }
}

# Validates cram files for formatting issues. 
task ValidateCram {
  input {
    File unmappedCram
    String base_file_name
    String docker
  }
  command {
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms2g" \
      ValidateSamFile \
        --INPUT ~{unmappedCram} \
        --MODE SUMMARY \
        --IGNORE_WARNINGS false > ~{base_file_name}.validation.txt
  }
  runtime {
    cpu: 2
    memory: "4 GB"
    docker: docker
  }
  output {
    File validation = "~{base_file_name}.validation.txt"
  }
}

task mergeBamstoCram {
  input {
    Array[File] bamsToMerge
    String base_file_name
    String docker
    Int threads
  }

  command {
    set -eo pipefail

    samtools merge -@ ~{threads-1} \
      --write-index --output-fmt CRAM  \
      ~{base_file_name}.merged.cram ~{sep=" " bamsToMerge} 

    }
  runtime {
    docker: docker
    cpu: threads
  }
  output {
    File cram = "~{base_file_name}.merged.cram"
    File crai = "~{base_file_name}.merged.cram.crai"
  }
}