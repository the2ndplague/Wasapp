#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process runPipeline {
    executor 'local'
    
    env RUN_NAME = params.run_name
    env WORKSPACE_DIR = params.workspace_dir
    env POD5_DIR = params.pod5_dir
    env DORADO_MODEL_PATH = params.dorado_model_path
    env KIT_NAME = params.kit_name
    env MIN_QSCORE = params.min_qscore
    env PRIMER_FILES_DIR = params.primer_files_dir
    env THREADS = params.threads
    env NOTEBOOKLM_EXEC = params.notebooklm_exec

    script:
    """
    cd "${params.workspace_dir}/scripts"
    chmod +x script_v2.sh
    ./script_v2.sh
    """
}

workflow {
    runPipeline()
}
