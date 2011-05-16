class Plotter {
	public:
		Plotter();
		void init(int argc, char **argv, int w, int h);
	private:
		static int loadData();
		
		static void renderScene(void);
		static void changeSize(int w, int h);
		static void mouseButton(int btn, int state, int x, int y);
		static void mouseMotion(int x, int y);
		static void specialKey(int key, int x, int y);

		static bool paintCylinder;
		static bool paintDisks;
		static int size;
		static float angle[2];
		static float translate[3];
		static int mouse_x, mouse_y, mouse_btn;
		static Plotter *instance;
};
