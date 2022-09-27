double CalcIntegral(double(*function)(double x), double a, double b, double steps_num) {
	double step = (b - a) / steps_num;
	double area = 0.0;
	for (int step_i = 0; step_i < steps_num; step_i++) {
		area += function(a + (step_i + 0.5) * step) * step;
	}
	return area;
}