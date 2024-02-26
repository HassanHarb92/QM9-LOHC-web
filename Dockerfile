FROM continuumio/miniconda3

# Copy your application code and environment.yml to the container
COPY . /app
WORKDIR /app

# Install your Conda environment
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "your_env_name", "/bin/bash", "-c"]

# The command to run your application
CMD conda run --no-capture-output -n your_env_name python app.py

