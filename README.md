# MetalMiner

## Dependencies

This project depends on the CSD Python API (`csd-python-api`), which is distributed via the CCDC conda channel.
Install it in a valid conda environment using a supported Python version:

```bash
conda install --channel=https://conda.ccdc.cam.ac.uk csd-python-api
```

You may also configure the channel once for future installs:

```bash
conda config --add channels https://conda.ccdc.cam.ac.uk
conda install csd-python-api
```

Note: a valid CCDC license is required for the API to function.

## Using in Jupyter

```python
from metalminer.main import main

# Preserve existing notebook workflow
result = main(
    TARGET_METAL_list=["Pu"],
    VISUALIZE=True,
    OXIDATION_STATES_FILTER=["all"],
)
```

You can also call the core API directly:

```python
from metalminer.core import Config, run_pipeline

config = Config(TARGET_METAL_list=["Pu"], VISUALIZE=True)
result = run_pipeline(config)
```

## Using the CLI

```bash
pip install -e .
python -m metalminer.cli --target-metals Pu --visualize --process-limit 100
```

## Running the Web GUI

```bash
pip install -e .
streamlit run app/app.py
```

The GUI mirrors all CLI options and streams console output during execution.
