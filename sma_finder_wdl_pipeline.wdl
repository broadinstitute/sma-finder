version 1.0

workflow SMAFinderWorkflow {
     input {
		 String genome_version = "hg38"   # can be "hg38", "hg37", or "t2t"
		 File ref_fasta
		 File ref_fasta_fai = "${ref_fasta}.fai"
		 File input_cram_or_bam
		 File input_crai_or_bai
		 String output_filename_prefix = sub(basename(input_cram_or_bam), "\\.bam$\|\\.cram$", "")
    }

	call SMAFinder {
        input:
			genome_version=genome_version,
            ref_fasta=ref_fasta,
            ref_fasta_fai=ref_fasta_fai,
			input_cram_or_bam=input_cram_or_bam,
			input_crai_or_bai=input_crai_or_bai,
			output_filename_prefix=output_filename_prefix
    }

	output {
		File output_tsv = SMAFinder.output_tsv
	}
}


task SMAFinder {

	input {
		String genome_version
    	File ref_fasta
		File ref_fasta_fai
		File input_cram_or_bam
		File input_crai_or_bai
		String output_filename_prefix

		Int disk_size = ceil(size(ref_fasta, "GB") + size(input_cram_or_bam, "GB") + 5)
	}

    command {
		if [[ "${genome_version}" != "hg37" && "${genome_version}" != "hg38" && "${genome_version}" != "t2t" ]]; then
			echo "ERROR: unexpected genome_version: ${genome_version}. It must be one of 'hg37', 'hg38', or 't2t'"
			exit 1
		fi

		set -ex
		python3 -u /sma_finder.py --verbose --${genome_version}-reference-fasta ${ref_fasta} \
			--output-tsv ${output_filename_prefix}.tsv ${input_cram_or_bam}
    }

    output {
		File output_tsv = "${output_filename_prefix}.tsv"
    }

     runtime {
         docker: "weisburd/sma_finder@sha256:41a3bc624a65d1625f098b7041ce6be310a7c002375ea43da3bc7e99cc9b00ed"
         cpu: 1
         preemptible: 1
         disks: "local-disk ${disk_size} HDD"
     }
}
