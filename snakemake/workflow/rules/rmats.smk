"""
rMATS Alternative Splicing Analysis
caveats: rMATS does not perform well as claimed, low statitical power for obvious AS events
"""

rule Read_Length_Detection:
    input:
        expand("bam_qc/stats/{sample}.stats.txt",sample=SAMPLES)
        if config['START'] in ["FASTQ", "BAM"]
        else "BAM file not provided. AND can't be generated from FASTQ, because FASTQ not provided"
    output:
        "meta/read_length.median.txt",
        "meta/read_length.max.txt"
    resources:
        mem_mb=1000,
    threads:
        1
    log:
        "meta/log/read_length.log"
    run:
        import re
        import statistics

        median_lens = []
        max_lens = []
        for f in input:
            p1 = re.compile("average length:\s*(\d*)")
            p2 = re.compile("maximum length:\s*(\d*)")
            for i, line in enumerate(open(f,"r")):
                for match1 in re.finditer(p1,line):
                    median_lens.append(match1.group(1))
                for match2 in re.finditer(p2,line):
                    max_lens.append(match2.group(1))
        median_lens = list(map(int,median_lens))
        max_lens = list(map(int,max_lens))
        print("average lengths detected from BAM:",median_lens,"\n")
        print("max lengths detected from BAM:",max_lens,"\n")
        if (len(set(max_lens)) < 1):
            sys.exit("read lengths detection from BAM failed")
        with open(output[0],"w") as out:
            out.write(str(int(statistics.median(median_lens))))
        with open(output[1],"w") as out:
            out.write(str(max(max_lens)))

if config['RMATS_ANALYSIS']:
    rule rMATS:
        #todo: update to docker: simple to install, faster, recommended on website
        input:
            b1=lambda wildcards: B1S[int(wildcards['ascn']) - 1],
            b2=lambda wildcards: B2S[int(wildcards['ascn']) - 1],
            gtf=config["GTF"],
            length_file="meta/read_length.median.txt",
            strand_file="meta/strandness.detected.txt",
        output:
            "rMATS.{ascn}/output/Results_JunctionCountsBased/SE.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/MXE.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A3SS.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A5SS.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/RI.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/SE.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/MXE.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A3SS.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A5SS.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/RI.MATS.JCEC.txt",
        conda:
            "../envs/rmats.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 10000,
        threads:
            4
        params:
            b1=lambda wildcards: B1[int(wildcards['ascn']) - 1],
            b2=lambda wildcards: B2[int(wildcards['ascn']) - 1],
            type="paired" if config['PAIR_END'] else "single",
            analysis="P" if config['PAIR_END'] else "U",
            length=get_read_length("meta/read_length.median.txt"),
            strandness=get_strandness("meta/strandness.detected.txt",config),
            MAX_FDR=config['MAX_FDR']
        log:
            "rMATS.{ascn}/rMATS.log"
        benchmark:
            "rMATS.{ascn}/rMATS.benchmark"
        shell:
            """
            rm -rf rMATS.{wildcards.ascn}
            mkdir -p rMATS.{wildcards.ascn}/
            mkdir -p rMATS.{wildcards.ascn}/params/
            echo {params.b1} > rMATS.{wildcards.ascn}/params/b1.txt
            echo {params.b2} > rMATS.{wildcards.ascn}/params/b2.txt

            rmats.py  \
            --b1 rMATS.{wildcards.ascn}/params/b1.txt \
            --b2  rMATS.{wildcards.ascn}/params/b2.txt \
            --gtf {input.gtf} \
            -t {params.type} \
            --readLength {params.length} \
            --variable-read-length \
            --libType {params.strandness} \
            --nthread {threads} \
            --tstat {threads} \
            --cstat 0.2 \
            --allow-clipping \
            --novelSS --mil 20 --mel 2000 \
            --od rMATS.{wildcards.ascn}/output/ \
            --tmp rMATS.{wildcards.ascn}/tmp/  > {log} 2>&1

            cp workflow/envs/rMATS_ReadMe.html rMATS.{wildcards.ascn}/output/  >>{log} 2>&1;
            rm -rf rMATS.{wildcards.ascn}/tmp/   >> {log} 2>&1;
            mkdir -p rMATS.{wildcards.ascn}/output/Results_JunctionCountsBased/ >> {log} 2>&1;
            mkdir -p rMATS.{wildcards.ascn}/output/Results_JunctionCountsAndExonCountsBased/  >> {log} 2>&1;
            mkdir -p rMATS.{wildcards.ascn}/output/fromGTF/ >> {log} 2>&1;
            mv rMATS.{wildcards.ascn}/output/*JCEC.txt rMATS.{wildcards.ascn}/output/Results_JunctionCountsAndExonCountsBased/ &>> {log}
            mv rMATS.{wildcards.ascn}/output/*JC.txt rMATS.{wildcards.ascn}/output/Results_JunctionCountsBased/ &>> {log}
            mv rMATS.{wildcards.ascn}/output/fromGTF*txt  rMATS.{wildcards.ascn}/output/fromGTF/ &>> {log}
            rm -f rMATS.{wildcards.ascn}/output/*raw.input* &>> {log}
            """

    rule rMATS_FLT:
        input:
            "rMATS.{ascn}/output/Results_JunctionCountsBased/SE.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/MXE.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A3SS.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A5SS.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/RI.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/SE.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/MXE.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A3SS.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A5SS.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/RI.MATS.JCEC.txt",
        output:
            "rMATS.{ascn}/output/Results_JunctionCountsBased/SE.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/MXE.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A3SS.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A5SS.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/RI.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/SE.MATS.JCEC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/MXE.MATS.JCEC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A3SS.MATS.JCEC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A5SS.MATS.JCEC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/RI.MATS.JCEC.sig.tsv",
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4000
        threads:
            1
        run:
            for f in input:
                df = pd.read_table(f)
                out = df[df['PValue'] < 1]
                outname = re.sub(".txt$",".sig.tsv",f)
                out.to_csv(outname,sep="\t")
            import shutil

            shutil.make_archive('rMATS.' + str(wildcards.ascn),'zip', \
                root_dir="./",base_dir='rMATS.' + str(wildcards.ascn))
            shutil.move('rMATS.' + str(wildcards.ascn) + '.zip','rMATS.' + str(wildcards.ascn))
