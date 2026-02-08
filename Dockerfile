FROM python:3.13-slim

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PIP_NO_CACHE_DIR=1

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt
RUN pip install --upgrade pip \
    && pip install -r /app/requirements.txt

COPY app /app/app
COPY templates /app/templates
COPY static /app/static
COPY README.md /app/README.md
COPY LICENSE /app/LICENSE

RUN mkdir -p /app/data/raw /app/data/processed

EXPOSE 8000

CMD ["gunicorn", "-w", "2", "-k", "gthread", "--threads", "4", "--timeout", "180", "-b", "0.0.0.0:8000", "app.main:create_app()"]
