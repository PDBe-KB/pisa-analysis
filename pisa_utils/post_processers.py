from abc import ABC, abstractmethod
import json
import logging

LOGGER = logging.getLogger(__name__)


class PostProcesser(ABC):
    def __init__(self, input_json_path, output_json_path):
        self.input_json_path = input_json_path
        self.output_json_path = output_json_path

        with open(self.input_json_path, "r") as f:
            self.data = json.load(f)

    @abstractmethod
    def parse(self):
        pass


class PostProcessComplexTable(PostProcesser):
    def parse(self):
        pass
