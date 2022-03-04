version 1.0

workflow cram {
  input {
    File? inputBam
    File? inputCram
    String outputFileNamePrefix = "output"
  }

  if (defined(inputBam) == true) {
    call bamToCram {
      input:
        inputBam = inputBam,
        outputFileNamePrefix = outputFileNamePrefix

    }
  }
  if (defined(inputCram) == true) {
    call cramToBam {
      input:
        inputCram = inputCram,
        outputFileNamePrefix = outputFileNamePrefix
    }
  }

   output {
        File? outputBam = cramToBam.outputBam
        File? outputBamIndex = cramToBam.outputBamIndex
        File? outputCram = bamToCram.outputCram
        File? outputCramIndex = bamToCram.outputCramIndex

    }


  parameter_meta {
    inputBam: "optional input bam file that will be converted to cram file"
    inputCram: "optional input cram file that will be converted to bam file"
    outputFileNamePrefix: "Prefix for output file"
  }

  meta {
    author: "Xuemei Luo"
    email: "xuemei.luo@oicr.on.ca"
    description: "Workflow for converting cram to bam or bam to cram"
    dependencies:
    [
      {
      name: "samtools/1.9",
      url: "https://github.com/samtools/samtools"
      }
    ]
  }
}


task bamToCram {
  input{
    File? inputBam
    String outputFileNamePrefix
    Int jobMemory = 12
    Int timeout = 5
    String referenceFasta = "$HG38_ROOT/hg38_random.fa"
    String? addParam
    String modules="samtools/1.9 hg38/p12"   
  }

  String resultCram = "~{outputFileNamePrefix}.cram"
  String resultCramIndex = "~{outputFileNamePrefix}.crai"

  command <<<
    set -euo pipefail

    #convert bam to cram
    samtools view -C -T ~{referenceFasta} ~{addParam} -o ~{resultCram} ~{inputBam}

    # index cram
    samtools index ~{resultCram} ~{resultCramIndex}


  >>>

  runtime {
    memory: "~{jobMemory}G"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File outputCram = "~{resultCram}"
    File outputCramIndex = "~{resultCramIndex}"
  }

  parameter_meta {
    inputBam: "Input bam file"
    outputFileNamePrefix: "Output prefix to prefix output file names with"
    jobMemory: "Memory (in GB) to allocate to the job"
    referenceFasta: "The fasta that is being used as a refrence to build the cram file"
    modules: "Modules required to process this step"
    timeout: "Hours before task timeout"
    addParam: "additional parameters"
  }

  meta {
    output_meta: {
      outputCram: "Converted cram file",
      outputCramIndex: "Index of the Converted CRAM file"
    }
  }
}

task cramToBam {
  input{
    File? inputCram
    String outputFileNamePrefix
    Int jobMemory = 12
    Int timeout = 5
    String referenceFasta = "$HG38_ROOT/hg38_random.fa"
    String? addParam
    String modules="samtools/1.9 hg38/p12"
  }

  String resultBam = "~{outputFileNamePrefix}.bam"
  String resultBamIndex = "~{outputFileNamePrefix}.bai"

  command <<<
    set -euo pipefail

    #convert cram to bam
    samtools view -b -T ~{referenceFasta} ~{addParam} -o ~{resultBam} ~{inputCram}

    # index bam
    samtools index ~{resultBam} ~{resultBamIndex}


  >>>

  runtime {
    memory: "~{jobMemory}G"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File outputBam = "~{resultBam}"
    File outputBamIndex = "~{resultBamIndex}"
  }

  parameter_meta {
    inputCram: "Input Cram file"
    outputFileNamePrefix: "Output prefix to prefix output file names with"
    jobMemory: "Memory (in GB) to allocate to the job"
    referenceFasta: "The fasta that is being used as a refrence to build the bam file"
    modules: "Modules required to process this step"
    timeout: "Hours before task timeout"
    addParam: "additional parameters"
  }

  meta {
    output_meta: {
      outputBam: "Converted BAM file",
      outputBamIndex: "Index of the Converted BAM file"
    }
  }
}
