# cramConverter

Workflow for converting cram to bam or bam to cram

## Overview

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)


## Usage

### Cromwell
```
java -jar cromwell.jar run cramConverter.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputAlign`|File|input alignment file that will be converted to cram or bam file
`isBam`|Boolean|is the input a bam


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|"output"|Prefix for output file


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`bamToCram.jobMemory`|Int|12|Memory (in GB) to allocate to the job
`bamToCram.timeout`|Int|5|Hours before task timeout
`bamToCram.referenceFasta`|String|"$HG38_ROOT/hg38_random.fa"|The fasta that is being used as a refrence to build the cram file
`bamToCram.addParam`|String?|None|additional parameters
`bamToCram.modules`|String|"samtools/1.9 hg38/p12"|Modules required to process this step
`cramToBam.jobMemory`|Int|12|Memory (in GB) to allocate to the job
`cramToBam.timeout`|Int|5|Hours before task timeout
`cramToBam.referenceFasta`|String|"$HG38_ROOT/hg38_random.fa"|The fasta that is being used as a refrence to build the bam file
`cramToBam.addParam`|String?|None|additional parameters
`cramToBam.modules`|String|"samtools/1.9 hg38/p12"|Modules required to process this step


### Outputs

Output | Type | Description
---|---|---
`outputAlign`|File|output alignment file
`outputIndex`|File|output alignment index file


## Commands
 This section lists command(s) run by cramConverter workflow
 
 * Running cram
 
 This workflow converts cram to bam or bam to cram
 
 <<<
     set -euo pipefail
 
     #convert bam to cram
     samtools view -C -T ~{referenceFasta} ~{addParam} -o ~{resultCram} ~{inputBam}
 
     # index cram
     samtools index ~{resultCram} ~{resultCramIndex}
 
 
   >>>
 <<<
     set -euo pipefail
 
     #convert cram to bam
     samtools view -b -T ~{referenceFasta} ~{addParam} -o ~{resultBam} ~{inputCram}
 
     # index bam
     samtools index ~{resultBam} ~{resultBamIndex}
 
 
   >>>
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
