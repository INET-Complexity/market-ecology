import esl

class SingleExperiment(esl.simulation.model):
    def __init__(self, environment, parameters):
        super().__init__(environment, parameters)
        self.economy = esl.simulation.world()


    def clear_market(initial_prices, excess_demand_functions):
        model = esl.economics.markets.walras.excess_demand_model(initial_prices)
        model.excess_demand_functions = excess_demand_functions
        return model.compute_clearing_quotes()

    def step(self, time_interval):
        """
            Here we do model-specific logic, such as deciding when to stop the simulation early
        :param time_interval:
        :return:
        """
        print(f"Here we do model-specific logic {time_interval}")
        return time_interval.upper



if __name__ == "__main__":
    environment = esl.computation.environment()
    parameters = esl.simulation.parameter.parametrization()

    model = SingleExperiment(environment, parameters)

    se = esl.simulation.time_interval(1,10)

    a = model.step(se)
    print(model.start)


