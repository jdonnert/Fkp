enum knot_types{
    KNOT_EMPTY,
    KNOT_START,          // = 1, 8bit x, 16bit y, 16bit dxdy, dxddy=0 
    KNOT_MINMAX,         // = 2, 3x8bit x, 16bit y0; y1=y2=0 
    KNOT_FULL,           // 3x8bit x, 3x16bit y
    KNOT_STOP,           // 8bit x, 16bit y, 16bit dxdy, dxddy=0
};

struct Knot {           // uncompressed knot
    enum knot_types type;  
    int8_t idx;         // knot index
    float P[2];         // knot coordinates
    float Mleft[2];     // left control point 
    float Mright[2];    // right control point 
	bool Is_Global_Max; 
};

#pragma pack(push) // push current packing val to stack
#pragma pack(1) // tell compiler not to use padding here
struct half_point { // 40 bit compressed, type HALF_KNOT
    uint8_t  x;
    uint16_t y;
    union control_points {
        uint8_t xLR[2];         // extrema have x values only 
		uint8_t xyLR[2];		// edges have only one control point
    } M;
};

struct full_point { // 72 bit compressed, type FULL_KNOT
    uint8_t x;
    uint16_t y;
    uint8_t xleft;
	uint8_t xright;
    uint16_t yleft;
    uint16_t yright;
};
#pragma pack(pop) // pop/restore former packing value

static void find_first_knots(const double *, struct Knot *),
            draw_curve(const struct Knot *, double *),
            update_control_points(struct Knot *);

static bool add_knot(const double *,struct Knot *, int, enum knot_types);

static float find_max_error(const struct Knot *,const double *,double *,int *);

static double compute_spec_parameters(double *);

static void compress_knots_binary(const struct Knot *, char *);
static int uncompress_knots_binary(const char*, struct Knot *);

uint16_t compressFloat_16bit(float); 
uint8_t compressFloat_8bit(float);  
float uncompressFloat_16bit(uint16_t), uncompressFloat_8bit(uint8_t);

