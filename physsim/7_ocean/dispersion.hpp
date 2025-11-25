namespace physsim
{
    /**
     * @brief Dispersion model that translates wave numbers "k" into wave speed c_p.
     */
    class Dispersion
    {
    public:
        /**
         * @brief Constructor.
         */
        Dispersion();

        /**
         * @brief Computes the phase speed c_p for a given wave number k
         * @param waveNumber wave number k
         * @return phase speed c_p [m/s]
         */
        double phaseSpeed(double waveNumber) const;

        /**
         * @brief Gravity [m/s^2]
         */
        double gravity;

        /**
         * @brief Water depth at the current point of evaluation. [m]
         */
        double depth;

        /**
         * @brief Liquid density rho [kg/m^3]
         */
        double density;

        /**
         * @brief Surface tension sigma [N/m]
         */
        double surfaceTension;
    };
}
