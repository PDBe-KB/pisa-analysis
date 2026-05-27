.PHONY: test

default: test

test:
	coverage run -m pytest --junitxml=report.xml -vv -s
	coverage xml -o coverage/cobertura-coverage.xml
	coverage report -m
