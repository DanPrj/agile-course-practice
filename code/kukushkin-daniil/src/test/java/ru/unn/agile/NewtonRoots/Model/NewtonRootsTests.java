package ru.unn.agile.NewtonRoots.Model;

import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.List;
import static java.lang.Double.NaN;

import org.junit.Test;

import static java.lang.Double.isNaN;
import static org.junit.Assert.*;

public class NewtonRootsTests {
    private NewtonMethod newtonMethod;
    private double intervalStart = 0, intervalEnd = 2, initialPoint = 1.5, eps = 1e-8, derivativeStep = 1e-3;

    @Test
    public void findRootIfExistInIntervalWithFunctionModuleStopping() {
        newtonMethod = new NewtonMethod(eps, derivativeStep);
        newtonMethod.setStoppingCriterion(NewtonMethod.StoppingCriterion.FunctionModule);
        FunctionInterface func = (x) -> (x - 1);

        double root = newtonMethod.findRoot(func, initialPoint, intervalStart, intervalEnd);

        assertTrue(Math.abs( func.compute(root) ) < eps);
    }

    @Test
    public void findRootIfExistInIntervalWithDifferenceBetweenApproxStopping() {
        newtonMethod = new NewtonMethod(eps, derivativeStep);
        newtonMethod.setStoppingCriterion(NewtonMethod.StoppingCriterion.DifferenceBetweenApproximates);
        FunctionInterface func = (x) -> (x - 1);

        double root = newtonMethod.findRoot(func, initialPoint, intervalStart, intervalEnd);

        assertTrue(Math.abs( func.compute(root) ) < eps);
    }

    @Test
    public void findRootOnTheLeftBorder() {
        newtonMethod = new NewtonMethod(eps, derivativeStep);
        FunctionInterface func = (x) -> (x);

        double root = newtonMethod.findRoot(func, initialPoint, intervalStart, intervalEnd);

        assertTrue(Math.abs( func.compute(root) ) < eps);
    }

    @Test
    public void findRootOnTheRightBorder() {
        newtonMethod = new NewtonMethod(eps, derivativeStep);
        FunctionInterface func = (x) -> (2 - x);

        double root = newtonMethod.findRoot(func, initialPoint, intervalStart, intervalEnd);

        assertTrue(Math.abs( func.compute(root) ) < eps);
    }

    @Test
    public void findRootIfNotExistInInterval() {
        newtonMethod = new NewtonMethod(eps, derivativeStep);
        FunctionInterface func = (x) -> (x + 1);

        double root = newtonMethod.findRoot(func, initialPoint, intervalStart, intervalEnd);

        assertTrue( isNaN(root) );
    }

    @Test
    public void findRootIfIncorrectIntervalBoundaries() {
        newtonMethod = new NewtonMethod(eps, derivativeStep);
        FunctionInterface func = (x) -> (x + 1);

        double root = newtonMethod.findRoot(func, initialPoint, intervalEnd, intervalStart);

        assertTrue( newtonMethod.resultStatus == NewtonMethod.ResultStatus.IncorrectIntervalBoundaries );
    }

    @Test
    public void findRootIfInitialPointOutsideInterval() {
        newtonMethod = new NewtonMethod(eps, derivativeStep);
        FunctionInterface func = (x) -> (x + 1);

        double root = newtonMethod.findRoot(func, -10, intervalStart, intervalEnd);

        assertTrue( newtonMethod.resultStatus == NewtonMethod.ResultStatus.InitialPointOutsideInterval );
    }

    @Test
    public void findRootIfNonmonotonicFunction() {
        newtonMethod = new NewtonMethod(eps, derivativeStep);
        FunctionInterface func = (x) -> (x - 1) * (x - 1);

        double root = newtonMethod.findRoot(func, initialPoint, intervalStart, intervalEnd);

        assertTrue( newtonMethod.resultStatus == NewtonMethod.ResultStatus.NonmonotonicFunctionOnInterval );
    }

    @Test
    public void setStoppingCriterionAsFunctionModule() {
        newtonMethod = new NewtonMethod();

        newtonMethod.setStoppingCriterion(NewtonMethod.StoppingCriterion.FunctionModule);

        assertTrue( newtonMethod.getStoppingCriterion() == NewtonMethod.StoppingCriterion.FunctionModule );
    }

    @Test
    public void setStoppingCriterionAsDifferenceBetweenApproximates() {
        newtonMethod = new NewtonMethod();

        newtonMethod.setStoppingCriterion(NewtonMethod.StoppingCriterion.DifferenceBetweenApproximates);

        assertTrue( newtonMethod.getStoppingCriterion() == NewtonMethod.StoppingCriterion.DifferenceBetweenApproximates );
    }

    @Test
    public void setCorrectAccuracyEps() {
        newtonMethod = new NewtonMethod();

        newtonMethod.setAccuracyEps(1e-4);

        assertTrue( newtonMethod.getAccuracyEps() == 1e-4);
    }

    @Test
    public void setZeroAccuracyEps() {
        newtonMethod = new NewtonMethod();

        newtonMethod.setAccuracyEps(0);

        assertTrue( newtonMethod.getAccuracyEps() == NewtonMethod.DEFAULT_EPS);
    }

    @Test
    public void setNegativeAccuracyEps() {
        newtonMethod = new NewtonMethod();

        newtonMethod.setAccuracyEps(-1);

        assertTrue( newtonMethod.getAccuracyEps() == NewtonMethod.DEFAULT_EPS);
    }


    @Test
    public void setCorrectDerivativeStep() {
        newtonMethod = new NewtonMethod();

        newtonMethod.setDerivativeStep(1e-2);

        assertTrue( newtonMethod.getDerivativeStep() == 1e-2);
    }

    @Test
    public void setZeroDerivativeStep() {
        newtonMethod = new NewtonMethod();

        newtonMethod.setDerivativeStep(0);

        assertTrue( newtonMethod.getDerivativeStep() == NewtonMethod.DEFAULT_DERIVATIVE_STEP);
    }

    @Test
    public void setNegativeDerivativeStep() {
        newtonMethod = new NewtonMethod();

        newtonMethod.setDerivativeStep(-1e-2);

        assertTrue( newtonMethod.getDerivativeStep() == NewtonMethod.DEFAULT_DERIVATIVE_STEP);
    }
}
