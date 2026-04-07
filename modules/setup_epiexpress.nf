nextflow.enable.dsl = 2

process SETUP_EPIEXPRESS {
    tag "EpiExpress"
    storeDir params.assets_dir ?: "assets"

    output:
    path "EpiExpress", emit: repo_dir

    script:
    def repo_url = 'https://github.com/SchulzLab/EpiExpress.git'
    def version = params.epiexpress_version ?: 'main'
    """
    set -e
    echo "=== EpiExpress Repository Setup ==="
    echo "Repository: ${repo_url}"
    echo "Version: ${version}"
    
    # Clone repository
    git clone --depth 1 --branch ${version} ${repo_url}

    if [ ! -d "EpiExpress" ]; then
    echo "ERROR: Clone failed" >&2
    exit 1
    fi
    """
}
