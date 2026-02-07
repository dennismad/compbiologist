PYTHON ?= python3

.PHONY: venv install run refresh test check

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
