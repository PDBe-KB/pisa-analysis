.PHONY: test

default: test

test:
	coverage run -m pytest --junitxml=report.xml -vv -s
	coverage xml -o coverage/cobertura-coverage.xml
	coverage report -m

lint:
	uv run pre-commit run --all

lest: lint test

requirements:
	uv export --no-dev --no-hashes --format requirements.txt -o requirements.txt

dev-requirements:
	uv export --only-dev --no-hashes --format requirements.txt -o requirements-dev.txt

all-requirements: requirements dev-requirements
	@echo "All requirements files generated."
