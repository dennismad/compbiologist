PYTHON ?= .venv/bin/python
ifeq ($(wildcard $(PYTHON)),)
PYTHON := python3
endif

.PHONY: venv install run refresh test check docker-up docker-down docker-logs

venv:
	$(PYTHON) -m venv .venv

install:
	. .venv/bin/activate && pip install -r requirements.txt

run:
	$(PYTHON) -m app.main

refresh:
	$(PYTHON) -m app.main --refresh

test:
	$(PYTHON) -m pytest

check:
	$(PYTHON) -m compileall app tests templates static

docker-up:
	docker compose up -d --build

docker-down:
	docker compose down

docker-logs:
	docker compose logs -f compbiologist
