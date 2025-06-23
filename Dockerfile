# Imagen base oficial con Python
FROM python:3.12-slim

# Evita prompts de configuración interactiva
ENV DEBIAN_FRONTEND=noninteractive

# Actualiza e instala herramientas de sistema y bioinformática
RUN apt-get update && apt-get install -y --no-install-recommends \
    clustalo \
    emboss \
    primer3 \
    build-essential \
    curl \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Crea un entorno virtual para instalar paquetes de Python
RUN python3 -m venv /venv

# Activa el entorno virtual e instala dependencias de Python
COPY requirements.txt .
RUN /venv/bin/pip install --upgrade pip
RUN /venv/bin/pip install -r requirements.txt

# Copia el resto del proyecto
COPY . /app
WORKDIR /app

# Asegura que se use el entorno virtual por defecto
ENV PATH="/venv/bin:$PATH"

# Comando por defecto al iniciar el contenedor
CMD ["bash"]
