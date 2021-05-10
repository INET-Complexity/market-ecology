import esl



class Fund(esl.economics.company):
    def __init__(self, identifier, jurisdiction = esl.law.jurisdictions.US):
        super().__init__(identifier, jurisdiction)
        pass


class NoiseTrader(Fund):

    def __init__(self, identifier):
        super().__init__(identifier)

    def act(self, esl.simulation.time_interval si, esl.simulation.seed):



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
        # Here we do model-specific logic

        print(f"step {time_interval}")
        return time_interval.upper



if __name__ == "__main__":
    nt = NoiseTrader(esl.simulation.identity([1,2]))



    environment = esl.computation.environment()
    parameters = esl.simulation.parameter.parametrization()

    model = SingleExperiment(environment, parameters)
    print(model)
    print(environment)




    se = esl.simulation.time_interval(1,10)

    print(  environment.run  )

    a = model.step(se)
    print(model.start)


