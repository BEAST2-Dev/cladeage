package beast.math.distributions;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class LogNormalImpl implements ContinuousDistribution {
    double m_fMean;
    double m_fStdDev;
    NormalDistributionImpl m_normal = new NormalDistributionImpl(0, 1);

    public LogNormalImpl(double fMean, double fStdDev) {
        setMeanAndStdDev(fMean, fStdDev);
    }
    public double getMean() {
    	return m_fMean;
    }
    public double getSigma() {
    	return m_fStdDev;
    }

    void setMeanAndStdDev(double fMean, double fStdDev) {
        m_fMean = fMean;
        m_fStdDev = fStdDev;
        m_normal.setMean(fMean);
        m_normal.setStandardDeviation(fStdDev);
    }

    @Override
    public double cumulativeProbability(double x) throws MathException {
        return m_normal.cumulativeProbability(Math.log(x));
    }

    @Override
    public double cumulativeProbability(double x0, double x1) throws MathException {
        return cumulativeProbability(x1) - cumulativeProbability(x0);
    }

    @Override
    public double inverseCumulativeProbability(double p) throws MathException {
        return Math.exp(m_normal.inverseCumulativeProbability(p));
    }

    @Override
    public double density(double fX) {
        return m_normal.density(Math.log(fX)) / fX;
    }

    @Override
    public double logDensity(double fX) {
        return m_normal.logDensity(Math.log(fX)) - Math.log(fX);
    }
} // class LogNormalImpl
