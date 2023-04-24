from functools import reduce


class LogicUtils:
    @staticmethod
    def cumulate_conditions(condition_list):
        if not condition_list:
            return None
        return reduce(lambda a, b: a & b, condition_list)