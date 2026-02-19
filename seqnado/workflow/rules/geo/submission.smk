from pathlib import Path
from seqnado.workflow.helpers.geo import get_files_for_symlink, get_symlinked_files


rule geo_samples_table:
    input:
        OUTPUT_DIR + "/protocol.txt",
    output:
        OUTPUT_DIR + "/geo_submission/samples_table.txt",
    params: 
        output=OUTPUT,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/samples_table.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/samples_table.tsv",
    message: "Generating samples table for GEO submission",
    run:
        from seqnado.outputs.files import GEOFiles
        source_files = [str(p) for p in OUTPUT.files if "/geo_submission/" not in str(p)]

        df = GEOFiles(make_geo_submission_files=True,
                      assay=OUTPUT.assay,
                      design=OUTPUT.design_dataframe,
                      sample_names=OUTPUT.sample_names,
                      config=OUTPUT.config,
                      processed_files=source_files
                      ).metadata

        df.to_csv(output[0], sep="\t", index=False)


rule geo_upload_instructions:
    input:
        OUTPUT_DIR + "/geo_submission/samples_table.txt",
    output:
        instructions=OUTPUT_DIR + "/geo_submission/upload_instructions.txt",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/geo_upload_instructions.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/geo_upload_instructions.tsv",
    message: "Creating GEO upload instructions",
    run:
        import importlib.resources
        import seqnado.data

        source = importlib.resources.files(seqnado.data) / 'geo_upload_instructions.txt'
        with open(source, 'r') as f:
            with open(output.instructions, 'w') as f_out:
                f_out.write(f.read())


rule geo_symlink:
    input:
        files=lambda wc: get_files_for_symlink(OUTPUT=OUTPUT, OUTPUT_DIR=OUTPUT_DIR, wc=wc),
    output:
        flag=temp(OUTPUT_DIR + "/geo_submission/.symlinks_created"),
    params:
        output=OUTPUT,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/geo_symlink.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/geo_symlink.tsv",
    message: "Creating symlinks for GEO submission files",
    run:
        from pathlib import Path
        from seqnado.outputs.files import GEOFiles

        # Exclude geo_submission files to avoid circular dependencies
        source_files = [str(p) for p in OUTPUT.files if "/geo_submission/" not in str(p)]

        geo_files = GEOFiles(make_geo_submission_files=True,
                             assay=OUTPUT.assay,
                             design=OUTPUT.design_dataframe,
                             sample_names=OUTPUT.sample_names,
                             config=OUTPUT.config,
                             processed_files=source_files,
                             output=OUTPUT)

        fastqs = geo_files.raw_files
        processed_files = geo_files.processed_data_files

        src = Path(OUTPUT_DIR + "/fastqs")
        dest = Path(OUTPUT_DIR + "/geo_submission")

        # Create symlinks for raw files
        for sample_name, fastq in fastqs.items():
            for fq in fastq:
                fq = Path(fq)
                src_file = (src / fq.name).absolute().resolve()
                dest_file = dest / fq.name
                if not dest_file.exists():
                    dest_file.symlink_to(src_file)

        # Create symlinks for processed files
        for _, row in processed_files.iterrows():
            src_file = Path(row['path']).absolute().resolve()
            dest_file = dest / row['output_file_name']
            if not dest_file.exists():
                dest_file.symlink_to(src_file)

        # Create flag file to indicate completion
        Path(output.flag).touch()


rule geo_md5sum:
    input:
        flag=OUTPUT_DIR + "/geo_submission/.symlinks_created",
    output:
        OUTPUT_DIR + "/geo_submission/md5sums.txt",
    params:
        geo_dir=OUTPUT_DIR + "/geo_submission",
        output_obj=OUTPUT,
        output_dir=OUTPUT_DIR,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/md5sum.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/md5sum.tsv",
    message: "Generating MD5 checksums for GEO submission files",
    run:
        from pathlib import Path

        # Get the list of files dynamically
        files = get_symlinked_files(OUTPUT=params.output_obj, OUTPUT_DIR=params.output_dir, wc=wildcards)

        # Get basenames
        basenames = [Path(f).name for f in files]

        # Change to geo directory and compute md5sums
        import subprocess
        result = subprocess.run(
            ["md5sum"] + basenames,
            cwd=params.geo_dir,
            capture_output=True,
            text=True,
            check=True
        )

        # Write output
        with open(output[0], 'w') as f:
            f.write(result.stdout)


rule geo_md5_table:
    input:
        md5sums=OUTPUT_DIR + "/geo_submission/md5sums.txt",
    output:
        raw=OUTPUT_DIR + "/geo_submission/raw_data_checksums.txt",
        processed=OUTPUT_DIR + "/geo_submission/processed_data_checksums.txt",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/geo_md5_table.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/geo_md5_table.tsv",
    message: "Creating separate MD5 checksum tables for raw and processed data",
    run:
        import pandas as pd 

        df = pd.read_csv(input.md5sums, sep=r"\s+", header=None, names=["md5sum", "file"])
        df["is_raw"] = df["file"].str.contains(".fastq.gz")

        df_raw = df[df["is_raw"]]
        df_processed = df[~df["is_raw"]]

        for outfile, df_sub in zip([output.raw, output.processed], [df_raw, df_processed]):
            df_sub.rename(columns={"file": "file name", "md5sum": "file checksum"}).to_csv(outfile, index=False, sep="\t")


rule geo_move_to_upload:
    input:
        flag=OUTPUT_DIR + "/geo_submission/.symlinks_created",
        validated=OUTPUT_DIR + "/geo_submission/.validated",
    output:
        outdir = directory(OUTPUT_DIR + f"/geo_submission/{ASSAY.clean_name}")
    params:
        output_obj=OUTPUT,
        output_dir=OUTPUT_DIR,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/move_to_upload.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/move_to_upload.tsv",
    message: "Moving files to final GEO upload directory",
    run:
        from pathlib import Path
        import shutil

        infiles = get_symlinked_files(
            OUTPUT=params.output_obj, 
            OUTPUT_DIR=params.output_dir, 
            wc=wildcards
        )
        Path(output.outdir).mkdir(parents=True, exist_ok=True)

        for f in infiles:
            shutil.copy2(f, output.outdir)


rule geo_remove_headers_for_security:
    input:
        flag=OUTPUT_DIR + "/geo_submission/.symlinks_created",
    output:
        validated=temp(OUTPUT_DIR + "/geo_submission/.validated"),
    params:
        output_obj=OUTPUT,
        output_dir=OUTPUT_DIR,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/geo/remove_headers.log",
    benchmark: OUTPUT_DIR + "/.benchmark/geo/remove_headers.tsv",
    message: "Removing headers from TSV files for GEO submission",
    run:
        from pathlib import Path

        # Get the list of files dynamically
        infiles = get_symlinked_files(
            OUTPUT=params.output_obj, 
            OUTPUT_DIR=params.output_dir, 
            wc=wildcards
        )

        for fn in infiles:
            path = Path(fn)
            if path.suffix == '.tsv':
                dest = path.with_suffix('.no_headers.tsv')
                
                with open(path, 'r') as f:
                    with open(dest, 'w') as f_out:
                        for line in f:
                            if line.startswith('#'):
                                continue
                            
                            f_out.write(line)
                path.unlink()
                dest.rename(path)
        
        Path(OUTPUT_DIR + "/geo_submission/.validated").touch()


localrules:
    geo_symlink,
    geo_md5_table,
    geo_samples_table,
    geo_upload_instructions,
    geo_move_to_upload,
    geo_remove_headers_for_security
