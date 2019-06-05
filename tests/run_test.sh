python3 /codebase/SCRNApipeline/generate_config.py --cellranger /codebase/cellranger-3.0.2/cellranger-cs/3.0.2/bin/ --markers /codebase/markers.yaml -sampleid test --build test
python3 /codebase/pipeline_basic.py --submit local
