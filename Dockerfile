FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    clustalo \
    emboss \
    primer3 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /venv

COPY requirements.txt .
RUN /venv/bin/pip install --upgrade pip
RUN /venv/bin/pip install -r requirements.txt

COPY . /app
WORKDIR /app

ENV PATH="/venv/bin:$PATH"

CMD ["bash"]
