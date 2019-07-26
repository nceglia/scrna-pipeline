python3 /codebase/SCRNApipeline/generate_config.py --reference /codebase/SCRNApipeline/tests/reference --data /codebase/SCRNApipeline/tests/data --sampleid test --cellranger /codebase/cellranger-3.0.2/cellranger-cs/3.0.2/bin/ --markers /codebase/markers.yaml --sampleid test --build test
python3 /codebase/SCRNApipeline/pipeline_basic.py --submit local
