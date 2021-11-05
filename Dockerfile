from pinellolab/stream

RUN pip install networkx pandas matplotlib scvelo numpy dash dash-bootstrap-components flask flask-cors qrcode scvr scanpy pillow gunicorn && apt update --allow-releaseinfo-change -y && apt install -y libnss3-tools wget && wget https://github.com/FiloSottile/mkcert/releases/download/v1.4.3/mkcert-v1.4.3-linux-amd64 && chmod +x mkcert-v1.4.3-linux-amd64 && ./mkcert-v1.4.3-linux-amd64 -install && ./mkcert-v1.4.3-linux-amd64 localhost

WORKDIR /app
ADD . /app

ENTRYPOINT ["/bin/bash", "/app/start.sh"] 
