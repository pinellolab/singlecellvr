
FROM continuumio/miniconda3

RUN conda install python=3.6 && \
    pip install scanpy scvr dash_bootstrap_components && \
    conda install -c conda-forge dash>=1.0.0 plotly>=3.10.0 \
                  gunicorn>=19.7.1 numpy>=1.14.0 pandas>=0.22.0 \
                  qrcode pillow && \
    conda install -c bioconda cufflinks 

WORKDIR /app

EXPOSE 8000

ENTRYPOINT [ "gunicorn", "--bind", "0.0.0.0:8000", "--pythonpath", "dash_app/apps/dash-singlecell-vr", "app:server", "--timeout", "300" ]
