

class Parameters():

    def __init__(self, parameters, config):
        self.parameters = parameters
        self.config = config
        self.calibration_dir = config["calibration_dir"]

    def get_solve_methods(self,object_type):
        method = self.parameters[object_type]["method"]
        linemeas_method = self.parameters[object_type]["linemeas_method"]
        methods = []
        if method:
            methods.append(method)
        if linemeas_method:
            methods.append(linemeas_method)
        return methods

    def get_solve_method(self, object_type):
        return self.parameters[object_type]["method"]

    def get_objects(self):
        return self.parameters["objects"]
