package ru.unn.agile.NewtonRoots.Model;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Double.NaN;

interface FunctionInterface {
    public double compute(double x);
}

public class NewtonMethod {
    public enum StoppingCriterion { FunctionModule, DifferenceBetweenApproximates }
    public enum ResultStatus { RootSuccessfullyFound, NoRootInInterval, NonmonotonicFunctionOnInterval,
                                InitialPointOutsideInterval, IncorrectIntervalBoundaries }
    private StoppingCriterion stoppingCriterion;
    public ResultStatus resultStatus;
    private double accuracyEps;
    private double derivativeStep;
    public static final double DEFAULT_EPS = 1e-10;
    public static final double DEFAULT_DERIVATIVE_STEP = 1e-10;
    private final StoppingCriterion DefaultStoppingCriterion = StoppingCriterion.FunctionModule;

    interface StoppingCriterionInterface {
        public boolean check(FunctionInterface func, double x, double xPrev, double accuracy);
    }
    StoppingCriterionInterface stoppingCriterionAsFunctionModule = (func, x, xPrev, eps) -> Math.abs(func.compute(x)) < eps;
    StoppingCriterionInterface stoppingCriterionAsDiffBetweenApprox = (func, x, xPrev, eps) -> Math.abs(x - xPrev) < eps;
    StoppingCriterionInterface currentStoppingCriterionFunction;

    public NewtonMethod()
    {
        setStoppingCriterion(DefaultStoppingCriterion);
        accuracyEps = DEFAULT_EPS;
        derivativeStep = DEFAULT_DERIVATIVE_STEP;
    }

    public NewtonMethod(double accuracy, double derivativeComputeStep)
    {
        setStoppingCriterion(DefaultStoppingCriterion);
        accuracyEps = accuracy;
        derivativeStep = derivativeComputeStep;
    }

    private boolean isMonotonicFunctionOnInterval(FunctionInterface func, double intervalStart, double intervalEnd)
    {
        double x = intervalStart;
        double xStep = 1e-5;
        double val = func.compute(x+xStep), dif;
        double oldVal = func.compute(x);
        double oldDif = val - oldVal;
        x += xStep;

        while(x <= intervalEnd)
        {
            x += xStep;
            val = func.compute(x);
            dif = val - oldVal;
            if(Math.signum(dif) != Math.signum(oldDif))
                return false;
            oldVal = val;
            oldDif = dif;
        }
        return true;
    }

    private boolean isMonotonicFunctionHasRoot(FunctionInterface func, double intervalStart, double intervalEnd)
    {
        return func.compute(intervalStart) * func.compute(intervalEnd) <= 0;
    }

    double findRoot(FunctionInterface func, double initialPoint, double intervalStart, double intervalEnd)
    {
        double x0 = initialPoint;
        double x = x0;
        double xPrev = x0;
        double h = derivativeStep;

        if(intervalEnd <= intervalStart)
        {
            resultStatus = ResultStatus.IncorrectIntervalBoundaries;
            return NaN;
        }

        if(initialPoint < intervalStart || initialPoint > intervalEnd)
        {
            resultStatus = ResultStatus.InitialPointOutsideInterval;
            return NaN;
        }

        if(isMonotonicFunctionOnInterval(func, intervalStart, intervalEnd) == false) {
            resultStatus = ResultStatus.NonmonotonicFunctionOnInterval;
            return NaN;
        }

        if(isMonotonicFunctionHasRoot(func, intervalStart, intervalEnd) == false) {
            resultStatus = ResultStatus.NoRootInInterval;
            return NaN;
        }

        do {
            x = x - func.compute(x) / (func.compute(x+h) - func.compute(x)) * h;
            while(x < intervalStart || x > intervalEnd)
                x = (x + xPrev) / 2;
            xPrev = x;
        } while ( !currentStoppingCriterionFunction.check(func, x, xPrev, accuracyEps) );
        resultStatus = ResultStatus.RootSuccessfullyFound;
        return x;
    }

    public void setStoppingCriterion(StoppingCriterion newStoppingCriterion)
    {
        stoppingCriterion = newStoppingCriterion;
        switch(stoppingCriterion) {
            case FunctionModule:
                currentStoppingCriterionFunction = stoppingCriterionAsFunctionModule;
                break;
            case DifferenceBetweenApproximates:
                currentStoppingCriterionFunction = stoppingCriterionAsDiffBetweenApprox;
                break;
        }
    }

    public StoppingCriterion getStoppingCriterion()
    {
        return stoppingCriterion;
    }

    public void setAccuracyEps(double accuracy)
    {
        if(accuracy > 0)
            accuracyEps = accuracy;
        else
            accuracyEps = DEFAULT_EPS;
    }

    public double getAccuracyEps()
    {
        return accuracyEps;
    }

    public void setDerivativeStep(double derivativeComputeStep)
    {
        if(derivativeComputeStep > 0)
            derivativeStep = derivativeComputeStep;
        else
            derivativeStep = DEFAULT_DERIVATIVE_STEP;
    }

    public double getDerivativeStep()
    {
        return derivativeStep;
    }


}
