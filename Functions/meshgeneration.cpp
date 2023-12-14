import std;

class MeshGeneration {
public:
    MeshGeneration() = default;
private:
    static constexpr double PI = std::numbers::pi;
    static constexpr double DOUBLEMAX = std::numeric_limits<double>::max();

    struct Edge {
        //Edge() = default;
        Edge(std::size_t aa, std::size_t bb) noexcept
            : a{ aa < bb ? aa : bb }, b{ aa < bb ? bb : aa } {}
        std::size_t  a{ 0 }, b{ 0 };
    };

    struct Para {
        //Para() = default;        
        Para(std::size_t ne, std::size_t c) noexcept
            : ne1{ ne }, c1{ c } {}

        bool isFull{ false };
        // position of triangle and point c
        std::size_t ne1{ 0 }, c1{ 0 };
        std::size_t ne2{ 0 }, c2{ 0 };

        //relative length
        double rl{ 0.0 };
    };

    struct EdgeHash {
        std::size_t operator()(const Edge& edge) const noexcept {
            return std::hash<std::size_t>{}(edge.a) ^ std::hash<std::size_t>{}(edge.b);
        }
    };

    struct EdgeEqual {
        bool operator()(const Edge& edge1, const Edge& edge2) const  noexcept {
            return (edge1.a == edge2.a && edge1.b == edge2.b);
        }
    };

    using EdgeParas = std::unordered_map<Edge, Para, EdgeHash, EdgeEqual>;
    using EdgeTable = std::unordered_set<Edge, EdgeHash, EdgeEqual>;

    struct Node {
        //Node() = default;	
        Node(std::size_t ii, double xx, double yy) noexcept
            : nn{ ii }, x{ xx }, y{ yy } {}
        std::size_t nn{ 0 };
        double x{ 0.0 }, y{ 0.0 }, dx{ 0.0 }, dy{ 0.0 };
        bool isInitial{ false }; bool isNarrow{ false };
        double r{ DOUBLEMAX };  double rmax{ DOUBLEMAX }; double rt{ 0.0 };

        // positions of boundaries the node belongs
        std::vector<std::size_t> nbs;

        // positions of nearby nodes 
        std::vector<std::size_t> nnbs;

        // positions of neighbouring nodes 
        std::vector<std::size_t> nnbis;

        // cavity edges around the node 
        std::vector<Edge> edges;
    };
    using Nodes = std::vector<Node>;

    struct FourNumber { std::size_t a{ 0 }, b{ 0 }, ne{ 0 }, c{ 0 }; };
    struct Element {
        //Element() = default;
        Element(std::size_t nene, std::size_t aa, std::size_t bb, std::size_t cc) noexcept
            : ne{ nene }, a{ aa }, b{ bb }, c{ cc },
            fournumber3{ FourNumber{a,b,ne,c}, FourNumber{b,c,ne,a}, FourNumber{c,a,ne,b} } { }
        std::size_t ne{ 0 }, a{ 0 }, b{ 0 }, c{ 0 };
        std::array<FourNumber, 3> fournumber3;
    };
    using Elements = std::vector<Element>;

    struct Boundary {
        //Boundary() = default;
        Boundary(std::size_t nbnb, std::size_t aa, std::size_t bb)
            : nb{ nbnb }, ns{ aa, bb } {}
        std::size_t nb{ 0 };  std::list<std::size_t> ns;
        double length{ 0.0 }; double rmax{ DOUBLEMAX };
    };
    using Boundaries = std::vector<Boundary>;

    using UnsignedDoubles = std::vector<std::pair<std::size_t, double>>;

    // Helper function to show an error message
    void error(const std::string& errormessage) {
        std::cerr << errormessage;
        throw std::runtime_error(errormessage);
    }

    void readGeometry(const std::filesystem::path& geofile, Nodes& nodes,
        Elements& elements, Boundaries& boundaries) {
        std::ifstream ist{ geofile, std::ios_base::binary };
        ist.exceptions(ist.exceptions() | std::ios_base::badbit);
        if (!ist) error("Open input file failed: " + geofile.string() + ".\n");

        std::uint64_t u64 = 0; void* pu64 = &u64;  char* pchu64 = static_cast<char*>(pu64);
        double f64 = 0; void* pf64 = &f64;  char* pchf64 = static_cast<char*>(pf64);

        ist.read(pchu64, 8); if (!ist) error("Read size of Nodes failed.\n");
        const std::size_t szns = u64;
        nodes.reserve(szns);
        for (std::size_t i = 0; i < szns; i++) {
            ist.read(pchf64, 8); if (!ist) error("Read x failed.\n");
            double x = f64;
            ist.read(pchf64, 8); if (!ist) error("Read y failed.\n");
            double y = f64;
            nodes.emplace_back(i, x, y);
        }

        ist.read(pchu64, 8); if (!ist) error("Read size of Elements failed.\n");
        const std::size_t szes = u64;
        elements.reserve(szes);
        for (std::size_t i = 0; i < szes; i++) {
            ist.read(pchu64, 8); if (!ist) error("Read a failed.\n");
            std::size_t a = u64;
            ist.read(pchu64, 8); if (!ist) error("Read b failed.\n");
            std::size_t b = u64;
            ist.read(pchu64, 8); if (!ist) error("Read c failed.\n");
            std::size_t c = u64;
            elements.emplace_back(i, a, b, c);
        }

        ist.read(pchu64, 8); if (!ist) error("Read size of Boundaries failed.\n");
        const std::size_t szbs = u64;
        boundaries.reserve(szbs);
        for (std::size_t i = 0; i < szbs; i++) {
            ist.read(pchu64, 8); if (!ist) error("Read a failed.\n");
            std::size_t a = u64;
            ist.read(pchu64, 8); if (!ist) error("Read b failed.\n");
            std::size_t b = u64;
            boundaries.emplace_back(i, a, b);
        }
    }

    double getThetaSumAB(std::size_t a, std::size_t b, std::size_t c1, std::size_t c2, const Nodes& nodes) {
        const Node& nodea{ nodes.at(a) };    const Node& nodeb{ nodes.at(b) };
        const Node& nodec1{ nodes.at(c1) }; const Node& nodec2{ nodes.at(c2) };

        const double xa{ nodea.x }, ya{ nodea.y }, xb{ nodeb.x }, yb{ nodeb.y };
        const double xc1{ nodec1.x }, yc1{ nodec1.y }, xc2{ nodec2.x }, yc2{ nodec2.y };
        const double xc1a{ xa - xc1 }, yc1a{ ya - yc1 }, xc1b{ xb - xc1 }, yc1b{ yb - yc1 },
            xc2a{ xa - xc2 }, yc2a{ ya - yc2 }, xc2b{ xb - xc2 }, yc2b{ yb - yc2 };
        const double theta1 = std::acos((xc1a * xc1b + yc1a * yc1b)
            / (std::hypot(xc1a, yc1a) * std::hypot(xc1b, yc1b)));
        const double theta2 = std::acos((xc2a * xc2b + yc2a * yc2b)
            / (std::hypot(xc2a, yc2a) * std::hypot(xc2b, yc2b)));
        return  theta1 + theta2;
    }

    void deleteThreeElementEdgesFromEdgeParas(EdgeParas& edgeparas, const Element& element) {
        for (auto [a, b, ne, c] : element.fournumber3) {
            const auto it = edgeparas.find(Edge{ a,b });
            if (it != edgeparas.end()) {
                Para& para = it->second;
                if (para.isFull) {
                    if (para.ne2 == ne && para.c2 == c)
                        para.isFull = false;
                    else if (para.ne1 == ne && para.c1 == c) {
                        para.isFull = false;
                        para.ne1 = para.ne2; para.c1 = para.c2;
                    }
                    else
                        error("Cannot find the ne and c from the para! \n");
                }
                else
                    edgeparas.erase(it);
            }
            else
                error("Cannot find the edge from edgeparas!\n");
        }
    }

    void addThreeElementEdgesToEdgeParas(EdgeParas& edgeparas, const Element& element) {
        for (const auto& [a, b, ne, c] : element.fournumber3) {
            const Edge edge{ a,b };
            const auto it = edgeparas.find(edge);
            if (it != edgeparas.end()) {
                Para& para = it->second;
                if (para.isFull)
                    error("Cannot add an edge anymore as the para is full!\n");
                else {
                    para.isFull = true;
                    para.ne2 = ne; para.c2 = c;
                }
            }
            else
                edgeparas.emplace(edge, Para{ ne,c });
        }
    }

    void updateMeshConnections(EdgeParas& edgeparas, Elements& elements, const Nodes& nodes) {
        using EdgeThetas = std::unordered_map<Edge, double, EdgeHash, EdgeEqual>;
        EdgeThetas edgethetas;
        for (const auto& [edge, para] : edgeparas) {
            if (para.isFull) {
                double theta = getThetaSumAB(edge.a, edge.b, para.c1, para.c2, nodes);
                if (theta > PI + 1e-10)
                    edgethetas[edge] = theta;
            }
        }

        //size_t i = 0;
        while (edgethetas.size() > 0) {
            auto itedgethetamax = std::max_element(std::execution::par_unseq,
                edgethetas.begin(), edgethetas.end(),
                [](const auto& edgetheta_a, const auto& edgetheta_b) {
                    return edgetheta_a.second < edgetheta_b.second; });

            const Edge edgemax{ itedgethetamax->first };
            const Para paramax = edgeparas.at(edgemax);
            const auto [a, b] {edgemax};
            const std::size_t ne1{ paramax.ne1 }; const std::size_t ne2{ paramax.ne2 };
            const std::size_t c1{ paramax.c1 };   const std::size_t c2{ paramax.c2 };
            edgethetas.erase(itedgethetamax);

            Element& element1 = elements.at(ne1);
            Element& element2 = elements.at(ne2);

            deleteThreeElementEdgesFromEdgeParas(edgeparas, element1);
            deleteThreeElementEdgesFromEdgeParas(edgeparas, element2);

            std::vector<Edge> edge0s;
            for (auto [aa, bb, nnee, cc] : element1.fournumber3) {
                Edge edge0{ aa,bb };
                if (edge0.a != a || edge0.b != b)
                    edge0s.emplace_back(std::move(edge0));
            }
            for (auto [aa, bb, nnee, cc] : element2.fournumber3) {
                Edge edge0{ aa,bb };
                if (edge0.a != a || edge0.b != b)
                    edge0s.emplace_back(std::move(edge0));
            }

            element1 = Element{ ne1, c1, c2, a };
            addThreeElementEdgesToEdgeParas(edgeparas, element1);

            element2 = Element{ ne2, c1, c2, b };
            addThreeElementEdgesToEdgeParas(edgeparas, element2);

            for (const auto& edge0 : edge0s) {
                const Para& para0 = edgeparas.at(edge0);
                if (para0.isFull) {
                    double theta0 = getThetaSumAB(edge0.a, edge0.b, para0.c1, para0.c2, nodes);
                    if (theta0 > PI + 1e-10)
                        edgethetas[edge0] = theta0;
                    else
                        edgethetas.erase(edge0);
                }
            }

            //if ((i++) % 20 == 0) {
            //    writeMeshRadii(nodes, elements, "meshradii.dat");
            //	continue;
            //}
        }
    }

    void readMeshFeature(const std::filesystem::path& meshfeaturefile,
        double& RMAX, double& GRAD, UnsignedDoubles& nnrs, UnsignedDoubles& nbrs) {
        std::ifstream ist{ meshfeaturefile };
        ist.exceptions(ist.exceptions() | std::ios_base::badbit);
        if (!ist) error("Open input file failed: " + meshfeaturefile.string() + ".\n");

        std::string word;

        ist >> word; if (!ist) error("Input word of Hmax failed.\n");
        ist >> RMAX; if (!ist) error("Input Hmax failed.\n");  RMAX /= 2.0;

        ist >> word; if (!ist) error("Input word of Hgrad failed.\n");
        ist >> GRAD; if (!ist) error("Input Hgrad failed.\n");

        std::size_t sz{ 0 }; std::size_t index{ 0 }; std::size_t n{ 0 }; double h{ 0.0 };

        ist >> word; if (!ist) error("Input word of Hvertex failed.\n");
        ist >> sz; if (!ist) error("Input size of Hvertex failed.\n"); nnrs.reserve(sz);
        for (std::size_t i = 0; i < sz; i++) {
            ist >> index; if (!ist) error("Input index of Hvertex failed.\n");
            ist >> n; if (!ist) error("Input n of Hvertex failed.\n");
            ist >> h; if (!ist) error("Input Hvertex failed.\n");
            nnrs.emplace_back(n, h / 2.0);
        }

        ist >> word; if (!ist) error("Input word of Hedge failed.\n");
        ist >> sz; if (!ist) error("Input size of Hvedge failed.\n"); nbrs.reserve(sz);
        for (std::size_t i = 0; i < sz; i++) {
            ist >> index; if (!ist) error("Input index of Hedge failed.\n");
            ist >> n; if (!ist) error("Input n of Hedge failed.\n");
            ist >> h; if (!ist) error("Input Hedge failed.\n");
            nbrs.emplace_back(n, h / 2.0);
        }
    }

    void initializeMeshParameters(Nodes& nodes, Boundaries& boundaries, EdgeParas& edgeparas,
        const Elements& elements, const double& RMAX, const UnsignedDoubles& nnrs, const UnsignedDoubles& nbrs) {
        for (auto& element : elements)
            addThreeElementEdgesToEdgeParas(edgeparas, element);

        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [&](Node& node) { node.isInitial = true; node.rmax = RMAX; });

        std::for_each(std::execution::par_unseq, boundaries.begin(), boundaries.end(),
            [&](Boundary& boundary) {  boundary.rmax = RMAX;
        const Node& nodea = nodes.at(boundary.ns.front());
        const Node& nodeb = nodes.at(boundary.ns.back());
        boundary.length = std::hypot(nodeb.x - nodea.x, nodeb.y - nodea.y);
            });

        for (const auto& [nn, rmax] : nnrs) {
            auto& node = nodes.at(nn);
            if (rmax < node.rmax)  node.rmax = rmax;
        }

        for (const auto& [nb, rmax] : nbrs) {
            Boundary& boundary = boundaries.at(nb);
            if (rmax < boundary.rmax) boundary.rmax = rmax;
            Node& nodea = nodes.at(boundary.ns.front());
            Node& nodeb = nodes.at(boundary.ns.back());
            if (rmax < nodea.rmax)  nodea.rmax = rmax;
            if (rmax < nodeb.rmax)  nodeb.rmax = rmax;
        }

        for (const auto& boundary : boundaries) {
            const std::size_t a = boundary.ns.front();
            const std::size_t b = boundary.ns.back();
            const std::size_t nb = boundary.nb;
            nodes.at(a).nbs.push_back(nb);
            nodes.at(b).nbs.push_back(nb);
        }

        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [](Node& node) { node.r = node.rmax; });
    }

    void writeMeshRadii(const Nodes& nodes, const Elements& elements,
        const std::filesystem::path& meshradiifile) {
        std::ofstream ost{ meshradiifile, std::ios_base::binary };
        ost.exceptions(ost.exceptions() | std::ios_base::badbit);
        if (!ost) error("Can't open input file " + meshradiifile.string() + ".\n");

        std::uint64_t u64 = 0; const void* pu64 = &u64;  const char* pchu64 = static_cast<const char*>(pu64);
        double f64 = 0; const void* pf64 = &f64;  const char* pchf64 = static_cast<const char*>(pf64);

        u64 = nodes.size();  ost.write(pchu64, 8);
        if (!ost) error("Write size of Nodes failed.\n");
        for (auto& node : nodes) {
            f64 = node.x; ost.write(pchf64, 8);
            if (!ost) error("Write x of node failed.\n");
            f64 = node.y; ost.write(pchf64, 8);
            if (!ost) error("Write y of node failed.\n");
            f64 = node.r; ost.write(pchf64, 8);
            if (!ost) error("Write r of node failed.\n");
        }

        u64 = elements.size();  ost.write(pchu64, 8);
        if (!ost) error("Write size of Elements failed.\n");
        for (auto& element : elements) {
            u64 = element.a; ost.write(pchu64, 8);
            if (!ost) error("Write a of element failed.\n");
            u64 = element.b; ost.write(pchu64, 8);
            if (!ost) error("Write b of element failed.\n");
            u64 = element.c; ost.write(pchu64, 8);
            if (!ost) error("Write c of element failed.\n");
        }
    }

    bool isInteriorEdgeAroundCorners(const Node& nodea, const Node& nodeb, const Boundaries& boundaries) {
        const std::vector<std::size_t>& nbas = nodea.nbs;
        const std::vector<std::size_t>& nbbs = nodeb.nbs;
        for (const auto nba : nbas) {
            const auto& nas = boundaries.at(nba).ns;
            for (const auto nbb : nbbs) {
                const auto& nbs = boundaries.at(nbb).ns;
                if (nas.front() != nbs.back() && nbs.front() != nas.back())
                    continue;
                else
                    return true;
            }
        }
        return false;
    }

    std::size_t getPositionOfSharedBoundary(const std::vector<std::size_t>& nbas,
        const std::vector<std::size_t>& nbbs) {
        const size_t sza = nbas.size(); const size_t szb = nbbs.size();
        if (sza > 2 || szb > 2)
            error("Wrong size for the nbas or nbbs!\n");

        for (const auto nba : nbas) {
            for (const auto nbb : nbbs)
                if (nba != nbb)
                    continue;
                else
                    return nba;
        }
        error("There is no shared boundary!\n");
        return 0;
    }

    std::vector<std::size_t> getPositionsOfThreeNeighbouringBoundaries(const std::size_t nbab,
        const Nodes& nodes, const Boundaries& boundaries) {
        auto& boundary = boundaries.at(nbab);

        auto& nodea = nodes.at(boundary.ns.front());
        auto& nodeb = nodes.at(boundary.ns.back());

        std::vector<std::size_t> nbabs; nbabs.reserve(3);
        for (auto nb : nodea.nbs) nbabs.push_back(nb);
        for (auto nb : nodeb.nbs) nbabs.push_back(nb);

        return nbabs;
    }

    bool isTriangularElementAroundCorners(const Node& nodea, const Node& nodeb, const Node& nodec,
        const Nodes& nodes, const Boundaries& boundaries) {
        const std::vector<std::size_t>& nbas = nodea.nbs;
        const std::vector<std::size_t>& nbbs = nodeb.nbs;
        const std::vector<std::size_t>& nbcs = nodec.nbs;

        const std::size_t nbab = getPositionOfSharedBoundary(nbas, nbbs);
        const auto nbabs = getPositionsOfThreeNeighbouringBoundaries(nbab, nodes, boundaries);

        for (const auto nbc : nbcs)
            for (const auto nbab0 : nbabs)
                if (nbc != nbab0)
                    continue;
                else
                    return true;
        return false;
    }

    double getAngle(double xa, double ya, double xb, double yb, double xi, double yi) noexcept {
        return std::acos(((xa - xi) * (xb - xi) + (ya - yi) * (yb - yi)) /
            (std::hypot(xa - xi, ya - yi) * std::hypot(xb - xi, yb - yi)));
    }

    double getHeightOfTriangle(double xa, double ya, double xb, double yb, double xi, double yi) noexcept {
        return std::abs((xa - xi) * (yb - yi) - (xb - xi) * (ya - yi))
            / std::hypot(xb - xa, yb - ya);
    }

    void updateRadiiOfBoundaryNodeBubblesAroundNarrowRegions(Nodes& nodes, const Boundaries& boundaries,
        const EdgeParas& edgeparas, const double GRAD) {
        for (const auto& [edge, para] : edgeparas) {
            const auto& [a, b] {edge};  Node& nodea{ nodes.at(a) }; Node& nodeb{ nodes.at(b) };
            const double& xa{ nodea.x }; const double& ya{ nodea.y }; double& ra{ nodea.r };
            const double& xb{ nodeb.x }; const double& yb{ nodeb.y }; double& rb{ nodeb.r };
            if (para.isFull) {
                if (!isInteriorEdgeAroundCorners(nodea, nodeb, boundaries)) {
                    const double rr = std::hypot(xb - xa, yb - ya) / 5.196;//  5.61 5.196
                    if (rr < ra)  ra = rr;     if (rr < rb)  rb = rr;
                }
            }
            else {
                Node& nodec{ nodes.at(para.c1) };
                if (!isTriangularElementAroundCorners(nodea, nodeb, nodec, nodes, boundaries)) {
                    const double& xc{ nodec.x }; const double& yc{ nodec.y }; double& rc{ nodec.r };
                    const double theta_a = getAngle(xc, yc, xb, yb, xa, ya);
                    const double theta_b = getAngle(xc, yc, xa, ya, xb, yb);
                    if (theta_a < (PI / 2.0 + 1e-8) && theta_b < (PI / 2.0 + 1e-8)) {
                        const double h = getHeightOfTriangle(xa, ya, xb, yb, xc, yc);
                        const double rrc = h / 5.196;
                        const double rra = rrc + (GRAD - 1.0) / (GRAD + 1.0)
                            * std::sqrt(std::pow(xa - xc, 2) + std::pow(ya - yc, 2) - std::pow(h, 2));
                        const double rrb = rrc + (GRAD - 1.0) / (GRAD + 1.0)
                            * std::sqrt(std::pow(xb - xc, 2) + std::pow(yb - yc, 2) - std::pow(h, 2));
                        if (rrc < rc)  rc = rrc;   if (rra < ra)  ra = rra;   if (rrb < rb)  rb = rrb;
                    }
                }
            }
        }
    }

    void labelBoundaryNodeBubblesAroundNarrowRegions(Nodes& nodes, const Boundaries& boundaries,
        const EdgeParas& edgeparas, const double GRAD) {
        // From here, isInitial (true) means boundary node bubble
        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [](Node& node) { node.isInitial = true; });

        //const double RATIO{ std::sqrt(GRAD) };
        const double RATIO{ 1.1 };
        for (const auto& [edge, para] : edgeparas) {
            const auto& [a, b] {edge};  Node& nodea{ nodes.at(a) }; Node& nodeb{ nodes.at(b) };
            const double xa{ nodea.x }; const double ya{ nodea.y }; const double ra{ nodea.r };
            const double xb{ nodeb.x }; const double yb{ nodeb.y }; const double rb{ nodeb.r };
            if (para.isFull) {
                if (!isInteriorEdgeAroundCorners(nodea, nodeb, boundaries)) {
                    const double rr = std::hypot(xb - xa, yb - ya) / 5.196;//  5.61 5.196                    
                    if (rr < RATIO * ra)  nodea.isNarrow = true; if (rr < RATIO * rb)  nodeb.isNarrow = true;
                }
            }
            else {
                Node& nodec{ nodes.at(para.c1) };
                if (!isTriangularElementAroundCorners(nodea, nodeb, nodec, nodes, boundaries)) {
                    const double xc{ nodec.x }; const double yc{ nodec.y }; const double rc{ nodec.r };
                    const double theta_a = getAngle(xc, yc, xb, yb, xa, ya);
                    const double theta_b = getAngle(xc, yc, xa, ya, xb, yb);
                    if (theta_a < (PI / 2.0 + 1e-8) && theta_b < (PI / 2.0 + 1e-8)) {
                        const double rrc = getHeightOfTriangle(xa, ya, xb, yb, xc, yc) / 5.196;
                        if (rrc < RATIO * rc)  nodec.isNarrow = true;
                    }
                }
            }
        }
    }

    bool isTwoLineSegmentsCrossed(const double xa, const double ya, const double xb, const double yb,
        const double xc, const double yc, const double xd, const double yd) noexcept {
        const double v1 = ((xc - xa) * (yb - ya) - (xb - xa) * (yc - ya)) *
            ((xd - xa) * (yb - ya) - (xb - xa) * (yd - ya));
        const double v2 = ((xa - xc) * (yd - yc) - (xd - xc) * (ya - yc)) *
            ((xb - xc) * (yd - yc) - (xd - xc) * (yb - yc));
        if (v1 > 0.0 || v2 > 0.0)
            return false;
        else
            return true;
    }

    bool isRequiredNewCavityElement(const Edge& edge, const std::size_t c, const std::size_t nn,
        const Nodes& nodes, const EdgeTable& edge_exterior_table) {
        const std::size_t a{ edge.a }; const std::size_t b{ edge.b };
        if (edge_exterior_table.contains(Edge{ a,c }) || edge_exterior_table.contains(Edge{ b,c }))
            return true;

        const Node& nodec = nodes.at(c); const double xc{ nodec.x }; const double yc{ nodec.y };
        const Node& nodenn = nodes.at(nn); const double xnn{ nodenn.x }; const double ynn{ nodenn.y };
        const Node& nodea = nodes.at(a); const double xa{ nodea.x }; const double ya{ nodea.y };
        const Node& nodeb = nodes.at(b); const double xb{ nodeb.x }; const double yb{ nodeb.y };
        if (isTwoLineSegmentsCrossed(xa, ya, xb, yb, xc, yc, xnn, ynn))
            return true;

        return false;
    }

    void setNumbersOfNearbyNodes(Nodes& nodes, const Elements& elements, const EdgeParas& edgeparas) {
        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [](Node& node) { node.nnbis.clear(); });
        for (const auto& [edge, para] : edgeparas) {
            auto& [a, b] {edge};
            nodes.at(a).nnbis.push_back(b);   nodes.at(b).nnbis.push_back(a);
        }

        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [&](Node& node) {
                const std::size_t nn{ node.nn };
                const std::vector<std::size_t>& nnbis = node.nnbis;

                std::unordered_set<std::size_t> ne_cavity_table;
                for (const std::size_t nnbi : nnbis) {
                    const Para& para = edgeparas.at(Edge{ nn, nnbi });
                    ne_cavity_table.insert(para.ne1);
                    if (para.isFull) {
                        ne_cavity_table.insert(para.ne2);
                    }
                }

                EdgeTable edge_exterior_table;
                for (const std::size_t ne : ne_cavity_table) {
                    for (const FourNumber& fournumber : elements.at(ne).fournumber3) {
                        const auto [itedge, isInserted] = edge_exterior_table.insert(Edge{ fournumber.a,fournumber.b });
                        if (!isInserted) {
                            edge_exterior_table.erase(itedge);
                        }
                    }
                }

                std::list<Edge> edge_exterior_list(std::from_range, edge_exterior_table);
                while (edge_exterior_list.size() > 0) {
                    auto itedge = edge_exterior_list.begin();
                    const Edge edge = *itedge;
                    edge_exterior_list.erase(itedge);
                    if (edge_exterior_table.contains(edge)) {
                        const Para& para = edgeparas.at(edge);
                        const std::size_t a{ edge.a }; const std::size_t b{ edge.b };
                        if (para.isFull) {
                            const bool isFirstElementContainedInCavity = ne_cavity_table.contains(para.ne1);
                            const std::size_t ne{ isFirstElementContainedInCavity ? para.ne2 : para.ne1 };
                            const std::size_t c{ isFirstElementContainedInCavity ? para.c2 : para.c1 };

                            if (isRequiredNewCavityElement(edge, c, nn, nodes, edge_exterior_table)) {
                                ne_cavity_table.insert(ne);
                                edge_exterior_table.erase(edge);
                                const Edge edge_ac{ a,c };
                                const auto [itedge_ac, isInserted_ac] = edge_exterior_table.insert(edge_ac);
                                if (isInserted_ac)
                                    edge_exterior_list.push_back(std::move(edge_ac));
                                else
                                    edge_exterior_table.erase(itedge_ac);

                                const Edge edge_bc{ b,c };
                                const auto [itedge_bc, isInserted_bc] = edge_exterior_table.insert(edge_bc);
                                if (isInserted_bc)
                                    edge_exterior_list.push_back(std::move(edge_bc));
                                else
                                    edge_exterior_table.erase(itedge_bc);
                            }
                        }
                    }
                }

                std::unordered_set<std::size_t> nnb_table;
                for (const auto ne : ne_cavity_table) {
                    const Element& element = elements.at(ne);
                    nnb_table.insert(element.a);
                    nnb_table.insert(element.b);
                    nnb_table.insert(element.c);
                }
                nnb_table.erase(nn);

                node.nnbs = std::vector<std::size_t>(std::from_range, nnb_table);
            });
    }        

    void updateRadiiOfBoundaryNodeBubblesAccordingToGradient(Nodes& nodes, const Elements& elements,
        const EdgeParas& edgeparas, const double GRAD) {
        setNumbersOfNearbyNodes(nodes, elements, edgeparas);

        std::list<std::reference_wrapper<Node>> noderefs;
        for (auto& node : nodes) noderefs.emplace_back(node);

        while (noderefs.size() > 1) {
            auto itnoderef = std::min_element(std::execution::par_unseq, noderefs.begin(), noderefs.end(),
                [&](const Node& nodea, const Node& nodeb) {
                    return nodea.r < nodeb.r;
                });

            const Node& nodea = *itnoderef;
            //delete min
            noderefs.erase(itnoderef);

            const double xa{ nodea.x }; const double ya{ nodea.y }; const double ra{ nodea.r };
            for (const std::size_t& b : nodea.nnbs) {
                Node& nodeb{ nodes.at(b) };  double& rb{ nodeb.r };
                if (rb > ra) {
                    const double xb{ nodeb.x }; const double yb{ nodeb.y };
                    const double rtemp{ ra + (GRAD - 1.0) / (GRAD + 1.0) * std::hypot(xb - xa, yb - ya) };
                    if (rtemp < rb) rb = rtemp;
                }
            }
        }
    }


    void insertNodeBubblesAtBoundaries(Boundaries& boundaries, Nodes& nodes,
        Elements& elements, EdgeParas& edgeparas, const double GRAD) {
        for (Boundary& boundary : boundaries) {
            auto& ns = boundary.ns;   const double length = boundary.length;
            const size_t nb = boundary.nb;

            double sumradii{ 0.0 };
            std::multimap<double, std::list<std::size_t>::iterator> ratioits;
            for (auto ita = ns.begin(), itb = std::next(ita); itb != ns.end(); ita++, itb++) {
                const Node& nodea{ nodes.at(*ita) }; const Node& nodeb{ nodes.at(*itb) };
                const double& xa{ nodea.x }; const double& ya{ nodea.y }; const double& ra{ nodea.r };
                const double& xb{ nodeb.x }; const double& yb{ nodeb.y }; const double& rb{ nodeb.r };
                auto ratio = (ra + rb) / std::hypot(xb - xa, yb - ya);
                ratioits.emplace(ratio, ita);
                sumradii += ra + rb;
            }
            double sumradiiT{ sumradii };
            double rmean{ 0.0 }; for (auto& n : ns) rmean += nodes.at(n).r; rmean /= ns.size();

            while (sumradii < length && sumradiiT < length) {
                const auto itratioitmin = ratioits.begin();
                auto ita{ itratioitmin->second }; auto  itb{ std::next(ita) };
                const std::size_t a{ *ita };  const std::size_t b{ *itb };
                const auto& para = edgeparas.at(Edge{ a, b });
                const std::size_t ne{ para.ne1 }; const std::size_t c{ para.c1 };
                auto& nodea{ nodes.at(a) }; auto& nodeb{ nodes.at(b) }; auto& nodec{ nodes.at(c) };
                const auto xa{ nodea.x }; const auto ya{ nodea.y }; const auto ra{ nodea.r };
                const auto xb{ nodeb.x }; const auto yb{ nodeb.y }; const auto rb{ nodeb.r };
                const auto xc{ nodec.x }; const auto yc{ nodec.y }; const auto rc{ nodec.r };
                double xi{ (xa + xb) / 2.0 }, yi{ (ya + yb) / 2.0 };
                const double rtemp1{ ra + (GRAD - 1.0) / (GRAD + 1.0) * std::hypot(xi - xa, yi - ya) };
                const double rtemp2{ rb + (GRAD - 1.0) / (GRAD + 1.0) * std::hypot(xi - xb, yi - yb) };
                const double rtemp3{ rc + (GRAD - 1.0) / (GRAD + 1.0) * std::hypot(xi - xc, yi - yc) };
                const double ri{ std::min({rtemp1, rtemp2, rtemp3, boundary.rmax}) };

                //ratio of gap for practical node bubbles
                const double rpp{ (length - sumradii) / sumradii };
                //ratio of gap for theoretical node bubbles
                const double rtt{ (length - sumradiiT) / sumradiiT };
                //ratio of gap/overlap for practical node bubbles after node bubble insertion
                const double rpi{ std::abs(length - sumradii - 2.0 * ri) / (sumradii + 2.0 * ri) };
                //ratio of gap/overlap for theoretical node bubbles after node bubble insertion
                const double rti{ std::abs(length - sumradiiT - 2.0 * rmean) / (sumradiiT + 2.0 * rmean) };
 
                if ((rpp > 0.25 && rtt > 0.25) || (rpp > rpi && rtt > rti)) {                    
                    ratioits.erase(itratioitmin);

                    Element& element_old = elements.at(ne);
                    deleteThreeElementEdgesFromEdgeParas(edgeparas, element_old);

                    std::size_t i{ nodes.size() };
                    Node& node_new = nodes.emplace_back(i, xi, yi);
                    node_new.r = ri; node_new.rmax = boundary.rmax;
                    node_new.nbs.emplace_back(nb);

                    element_old = Element{ ne, a, i , c };
                    addThreeElementEdgesToEdgeParas(edgeparas, element_old);

                    const Element& element_new = elements.emplace_back(elements.size(), i, b, c);
                    addThreeElementEdgesToEdgeParas(edgeparas, element_new);

                    auto  iti = ns.insert(itb, i);
                    double ratioa = (ra + ri) / std::hypot(xi - xa, yi - ya);
                    ratioits.emplace(ratioa, ita);
                    double ratioi = (rb + ri) / std::hypot(xb - xi, yb - yi);
                    ratioits.emplace(ratioi, iti);

                    sumradii += 2.0 * ri;
                    sumradiiT += 2.0 * rmean;
                    //writeMeshRadii(nodes, elements, "meshradii.dat");
                    //continue;
                }
                else
                    break;
            }
            //writeMeshRadii(nodes, elements, "meshradii.dat");
        }
    }

    void updatePositionsOfInsertedBoundaryNodeBubbles(Nodes& nodes, const Boundaries& boundaries) {
        for (const Boundary& boundary : boundaries) {
            auto& ns = boundary.ns;
            if (ns.size() != 2) {
                double length = boundary.length;
                const std::size_t a{ ns.front() };  const std::size_t b{ ns.back() };
                const Node& nodea{ nodes.at(a) }; const Node& nodeb{ nodes.at(b) };
                const double& xa{ nodea.x }; const double& ya{ nodea.y }; const double& ra{ nodea.r };
                const double& xb{ nodeb.x }; const double& yb{ nodeb.y }; const double& rb{ nodeb.r };

                double sumr{ -(ra + rb) };  for (const std::size_t& n : ns)	sumr += 2.0 * nodes.at(n).r;
                const double scale = length / sumr;

                double lengthai = scale * ra;
                auto itbegin = std::next(ns.begin()), itend = std::prev(ns.end());
                for (auto iti = itbegin; iti != itend; iti++) {
                    const std::size_t  i{ *iti }; Node& nodei{ nodes.at(i) };
                    double& xi{ nodei.x }; double& yi{ nodei.y }; const double ri{ nodei.r };
                    lengthai += scale * ri; double lengthbi{ length - lengthai };
                    xi = (lengthbi * xa + lengthai * xb) / length;
                    yi = (lengthbi * ya + lengthai * yb) / length;
                    lengthai += scale * ri;
                }
            }
        }
    }

    void resetRadiiOfBoundaryNodeBubblesAccordingToBoundaryEdgeLength(Nodes& nodes, const Elements& elements,
        const EdgeParas& edgeparas, const double GRAD) {
        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [](Node& node) { node.rt = 0.0; });

        for (const auto& [edge, para] : edgeparas) {
            if (!para.isFull) {
                Node& nodea{ nodes.at(edge.a) }; Node& nodeb{ nodes.at(edge.b) };
                const double radius = std::hypot(nodeb.x - nodea.x, nodeb.y - nodea.y) / 4.0;
                nodea.rt += radius; nodeb.rt += radius;
            }
        }

        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [](Node& node) { if (node.rt > node.rmax) node.r = node.rmax;
            else if (node.rt > 0.0) node.r = node.rt; });
    }

    void updateRadiiOfBoundaryNodeBubbles(Nodes& nodes, const Elements& elements,
        const Boundaries& boundaries, const EdgeParas& edgeparas, const double GRAD) {
        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [](Node& node) { if (!node.isInitial) node.r = node.rmax; });
        updateRadiiOfBoundaryNodeBubblesAroundNarrowRegions(nodes, boundaries, edgeparas, GRAD);
        //writeMeshRadii(nodes, elements, "meshradii.dat");
        updateRadiiOfBoundaryNodeBubblesAccordingToGradient(nodes, elements, edgeparas, GRAD);
        //writeMeshRadii(nodes, elements, "meshradii.dat");
    }

    void setNumbersOfCavityElements(std::vector<std::size_t>& nes, const double xi, const double yi,
        const Nodes& nodes, const Elements& elements, const EdgeParas& edgeparas) {
        std::unordered_set<std::size_t> netable(std::from_range, nes);
        for (std::size_t ii = 0; ii < nes.size(); ii++) {
            for (const auto& [aa1, bb1, nee1, cc1] : elements.at(nes.at(ii)).fournumber3) {
                const Edge edge{ aa1, bb1 }; const Para& para = edgeparas.at(edge);
                if (para.isFull) {
                    const std::size_t nee2 = (para.ne1 != nee1) ? para.ne1 : para.ne2;
                    const auto itnee2 = netable.find(nee2);
                    if (itnee2 != netable.end())
                        continue;
                    else {
                        const Element& element2 = elements.at(nee2);
                        const Node& nodea{ nodes.at(element2.a) };
                        const double xa{ nodea.x }; const double ya{ nodea.y };
                        const Node& nodeb{ nodes.at(element2.b) };
                        const double xb{ nodeb.x }; const double yb{ nodeb.y };
                        const Node& nodec{ nodes.at(element2.c) };
                        const double xc{ nodec.x }; const double yc{ nodec.y };

                        // circumcenter and radius of Circumcircle
                        const double xcc = (((pow(xb, 2) - pow(xa, 2)) + (pow(yb, 2) - pow(ya, 2))) * (yc - ya)
                            - ((pow(xc, 2) - pow(xa, 2)) + (pow(yc, 2) - pow(ya, 2))) * (yb - ya))
                            / (2 * ((xb - xa) * (yc - ya) - (xc - xa) * (yb - ya)));

                        const double ycc = (((pow(xb, 2) - pow(xa, 2)) + (pow(yb, 2) - pow(ya, 2))) * (xc - xa)
                            - ((pow(xc, 2) - pow(xa, 2)) + (pow(yc, 2) - pow(ya, 2))) * (xb - xa))
                            / (2 * ((xc - xa) * (yb - ya) - (xb - xa) * (yc - ya)));

                        const double rcc = hypot(xcc - xa, ycc - ya);

                        if (hypot(xi - xcc, yi - ycc) < rcc) {
                            netable.insert(nee2);   nes.push_back(nee2);
                        }
                    }
                }
            }
        }
    }

    EdgeTable getCavityExteriorEdges(const std::vector<std::size_t> nes,
        const Elements& elements, const EdgeParas& edgeparas) {
        EdgeTable edgetable;
        for (const std::size_t ne : nes) {
            for (const auto& fournumber : elements.at(ne).fournumber3) {
                const auto [itedge, isInserted] = edgetable.insert(Edge{ fournumber.a, fournumber.b });
                if (!isInserted)
                    edgetable.erase(itedge);
            }
        }
        return edgetable;
    }

    void insertNodeAtMiddleEdge(Nodes& nodes, Elements& elements, EdgeParas& edgeparas,
        const Edge& edgemax, const Para& paramax, const double GRAD, const double RMAX) {

        const auto [a, b] {edgemax};
        const std::size_t ne1{ paramax.ne1 }; const std::size_t ne2{ paramax.ne2 };
        //const std::size_t c1{ paramax.c1 };   const std::size_t c2{ paramax.c2 };

        const Node& nodea{ nodes.at(a) }; const Node& nodeb{ nodes.at(b) };
        const double xi = (nodea.x + nodeb.x) / 2.0;
        const double yi = (nodea.y + nodeb.y) / 2.0;

        const std::size_t i{ nodes.size() };
        Node& nodei = nodes.emplace_back(i, xi, yi); double& ri = nodei.r; ri = RMAX;

        std::vector<std::size_t> nes{ ne1, ne2 };
        setNumbersOfCavityElements(nes, xi, yi, nodes, elements, edgeparas);
        const EdgeTable edgetable = getCavityExteriorEdges(nes, elements, edgeparas);

        for (const std::size_t ne : nes)
            deleteThreeElementEdgesFromEdgeParas(edgeparas, elements.at(ne));

        auto itne = nes.begin();
        for (const Edge& edge : edgetable) {
            if (itne != nes.end()) {
                const std::size_t ne = *itne;
                Element& element = elements.at(ne);
                element = Element{ ne, edge.a, edge.b, i };
                addThreeElementEdgesToEdgeParas(edgeparas, element);
                itne++;
            }
            else {
                const std::size_t ne = elements.size();
                const Element& element = elements.emplace_back(ne, edge.a, edge.b, i);
                addThreeElementEdgesToEdgeParas(edgeparas, element);
            }
        }

        std::unordered_set<std::size_t> js;
        for (const Edge& edge : edgetable) {
            js.insert(edge.a); js.insert(edge.b);
        }

        for (const std::size_t j : js) {
            const Node& nodej{ nodes.at(j) };
            const double rj = nodej.r;
            if (rj < ri) {
                const double xj{ nodej.x }; const double yj{ nodej.y };
                double rtemp{ rj + (GRAD - 1.0) / (GRAD + 1.0) * std::hypot(xi - xj, yi - yj) };
                if (rtemp < ri) ri = rtemp;
            }
        }

        std::for_each(std::execution::par_unseq, js.begin(), js.end(),
            [&](const std::size_t j) { const Edge edge{ i, j }; auto& para = edgeparas.at(edge);
        const Node& nodea{ nodes.at(edge.a) }; const Node& nodeb{ nodes.at(edge.b) };
        para.rl = std::hypot(nodeb.x - nodea.x, nodeb.y - nodea.y) / (nodea.r + nodeb.r);
            });

        //writeMeshRadii(nodes, elements, "meshradii.dat");
    }

    bool isNewPositionInCavity(const Node& node, const Nodes& nodes) {
        const double x = node.x; const double y = node.y;
        const double xn = x + node.dx; const double yn = y + node.dy;

        for (const auto& [a, b] : node.edges) {
            const Node& nodea{ nodes.at(a) }; const Node& nodeb{ nodes.at(b) };
            const double xa{ nodea.x }, ya{ nodea.y }, xb{ nodeb.x }, yb{ nodeb.y };

            const double area1 = (yb - y) * (xa - x) - (xb - x) * (ya - y);
            const double area2 = (yb - yn) * (xa - xn) - (xb - xn) * (ya - yn);

            if (area1 * area2 < 0.0)
                return false;
        }

        return true;
    }

    void updatePositionsOfInsertedInteriorNodeBubbles(Nodes& nodes, const EdgeParas& edgeparas, const double ALPHA) {
        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [](Node& node) { node.nnbis.clear(); node.edges.clear(); });
        for (const auto& [edge, para] : edgeparas) {
            const auto& [a, b] {edge};
            nodes.at(a).nnbis.push_back(b);
            nodes.at(b).nnbis.push_back(a);
            if (para.isFull) {
                nodes.at(para.c1).edges.push_back(edge);
                nodes.at(para.c2).edges.push_back(edge);
            }
            else
                nodes.at(para.c1).edges.push_back(edge);
        }

        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [&](Node& node) {                
                if (!node.isInitial) {
                    double& dx = node.dx; double& dy = node.dy;
                    dx = 0.0; dy = 0.0;
                    for (const std::size_t a : node.nnbis) {
                        const Node& nodea{ nodes.at(a) };
                        const double length = std::hypot(nodea.x - node.x, nodea.y - node.y);
                        const double r = ALPHA * (node.r + nodea.r);
                        const double f = (length > r) ? 0.0 : (length - r);
                        const double nx = (nodea.x - node.x) / length;
                        const double ny = (nodea.y - node.y) / length;
                        dx += 0.3 * f * nx;  dy += 0.3 * f * ny;
                    }
                    
                    while (!isNewPositionInCavity(node, nodes) ) {
                        dx /= 2.0; dy /= 2.0;            
                    }
                }
            });

        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [](Node& node) { if(!node.isInitial) {node.x += node.dx; node.y += node.dy; }});
    }

    void updateRadiiOfInteriorNodeBubbles(Nodes& nodes, const Elements& elements,
        const EdgeParas& edgeparas, const double GRAD, const double RMAX) {
        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [](Node& node) { node.nnbis.clear(); });
        for (const auto& [edge, para] : edgeparas) {
            auto& [a, b] {edge};
            nodes.at(a).nnbis.push_back(b);   nodes.at(b).nnbis.push_back(a);
        }

        std::for_each(std::execution::par_unseq, nodes.begin(), nodes.end(),
            [&](Node& node) {
                if (!node.isInitial) {
                    bool isAroundNarrow = false; 
                    double sumradii{ 0.0 }; std::size_t szsum{ 0 };
                    for (const auto nnbi : node.nnbis) {
                        const Node& nodenbi = nodes.at(nnbi); 
                        if (nodenbi.isInitial) { sumradii += nodenbi.r; ++szsum; }
                        if (nodenbi.isNarrow) isAroundNarrow = true;
                        //if (nodenbi.isNarrow) { sumradii += nodenbi.r; ++szsum; }
                    }
                    //if (szsum > 1) isAroundNarrow = true;
                    node.r = isAroundNarrow ? (sumradii / szsum) : RMAX;
                }
            });

        std::list<std::reference_wrapper<Node>> noderefs;
        for (auto& node : nodes) noderefs.emplace_back(node);
        while (noderefs.size() > 1) {
            auto itnoderef = std::min_element(std::execution::par_unseq, noderefs.begin(), noderefs.end(),
                [](const Node& nodea, const Node& nodeb)
                { return nodea.r < nodeb.r; });

            const Node& nodea = *itnoderef;
            //delete min
            noderefs.erase(itnoderef);

            const double xa{ nodea.x }; const double ya{ nodea.y }; const double ra{ nodea.r };
            for (const std::size_t b : nodea.nnbis) {
                auto& nodeb{ nodes.at(b) };  auto& rb{ nodeb.r };
                if (!nodeb.isInitial && rb > ra) {
                    const double xb{ nodeb.x }; const double yb{ nodeb.y };
                    double rtemp{ ra + (GRAD - 1.0) / (GRAD + 1.0) * std::hypot(xb - xa, yb - ya) };
                    if (rtemp < rb) rb = rtemp;
                }
            }
        }
    }

    void updateRelativeLengthOfEdges(const Nodes& nodes, EdgeParas& edgeparas) noexcept {
        std::for_each(std::execution::par_unseq, edgeparas.begin(), edgeparas.end(),
            [&](auto& edgepara) { auto& [edge, para] = edgepara;
        if (para.isFull) {
            const Node& nodea{ nodes.at(edge.a) }; const Node& nodeb{ nodes.at(edge.b) };
            para.rl = std::hypot(nodeb.x - nodea.x, nodeb.y - nodea.y) / (nodea.r + nodeb.r);
        }
            });
    }

    void generateInteriorMeshes(Nodes& nodes, Elements& elements, EdgeParas& edgeparas,
        double GRAD, double RMAX) {
        updateRelativeLengthOfEdges(nodes, edgeparas);

        auto start = std::chrono::steady_clock::now();
        for (std::size_t ii = 0; ; ii++) {
            auto itedgeparamax = std::max_element(std::execution::par_unseq,
                edgeparas.begin(), edgeparas.end(),
                [](const auto& edgepara_a, const auto& edgepara_b)
                { return edgepara_a.second.rl < edgepara_b.second.rl; });
            const Edge& edgemax = itedgeparamax->first;   const Para& paramax = itedgeparamax->second;
            const double rlmax = paramax.rl;

            if (rlmax > 1.4 && paramax.isFull) {
                insertNodeAtMiddleEdge(nodes, elements, edgeparas, edgemax, paramax, GRAD, RMAX);
                //writeMeshRadii(nodes, elements, "meshradii.dat"); 
                if (rlmax < 1.6) {
                    //writeMeshRadii(nodes, elements, "meshradii.dat");
                    for (std::size_t ii = 0; ii < 3; ii++) {                    
                        for (std::size_t i = 0; i < 10; i++) {
                            updatePositionsOfInsertedInteriorNodeBubbles(nodes, edgeparas, 1.0);
                            //writeMeshRadii(nodes, elements, "meshradii.dat");

                            updateMeshConnections(edgeparas, elements, nodes);
                            //writeMeshRadii(nodes, elements, "meshradii.dat"); 
                        }
                        updateRadiiOfInteriorNodeBubbles(nodes, elements, edgeparas, GRAD, RMAX);
                        //writeMeshRadii(nodes, elements, "meshradii.dat");
                    }      

                    updateRelativeLengthOfEdges(nodes, edgeparas);
                    //writeMeshRadii(nodes, elements, "meshradii.dat");                    
                }

                const auto end = std::chrono::steady_clock::now();
                const std::chrono::duration<double> elapsed_seconds = end - start;
                if (elapsed_seconds.count() > 15.0) {
                    std::cout << "Iteration in progress, number of current nodes: " << nodes.size() << '\n';
                    start = std::chrono::steady_clock::now();
                }
            }
            else
                break;
        }

        //writeMeshRadii(nodes, elements, "meshradii.dat");

        for (std::size_t ii = 0; ii < 30; ii++) {
            for (std::size_t i = 0; i < 10; i++) {
                updatePositionsOfInsertedInteriorNodeBubbles(nodes, edgeparas, 1.2);
                //writeMeshRadii(nodes, elements, "meshradii.dat");

                updateMeshConnections(edgeparas, elements, nodes);
                //writeMeshRadii(nodes, elements, "meshradii.dat");
            }

            //writeMeshRadii(nodes, elements, "meshradii.dat");
            updateRadiiOfInteriorNodeBubbles(nodes, elements, edgeparas, GRAD, RMAX);
            //writeMeshRadii(nodes, elements, "meshradii.dat");
        }
    }

public:
    void implement(const std::filesystem::path& geofile, const std::filesystem::path& meshfeaturefile,
        const std::filesystem::path& meshradiifile) {

        //Step 0: Initialization of Bubbles and Predefined Parameters

        Nodes nodes; Elements elements; Boundaries boundaries; EdgeParas edgeparas;
        readGeometry(geofile, nodes, elements, boundaries);
        //writeMeshRadii(nodes, elements, "meshradii.dat");

        updateMeshConnections(edgeparas, elements, nodes);
        //writeMeshRadii(nodes, elements, "meshradii.dat");

        double RMAX{ 0.0 }; double GRAD{ 0.0 };
        UnsignedDoubles nnrs;  UnsignedDoubles nbrs;

        readMeshFeature(meshfeaturefile, RMAX, GRAD, nnrs, nbrs);

        initializeMeshParameters(nodes, boundaries, edgeparas, elements, RMAX, nnrs, nbrs);
        //writeMeshRadii(nodes, elements, "meshradii.dat");

        //Step 1: Refinement of Mesh Boundaries

        std::size_t round{ 0 };
        while (true) {
            const std::size_t sz_old{ nodes.size() };
                                     
            resetRadiiOfBoundaryNodeBubblesAccordingToBoundaryEdgeLength(nodes, elements, edgeparas, GRAD);
            //writeMeshRadii(nodes, elements, "meshradii.dat");

            updateRadiiOfBoundaryNodeBubblesAccordingToGradient(nodes, elements, edgeparas, GRAD);
            //writeMeshRadii(nodes, elements, "meshradii.dat");

            updateRadiiOfBoundaryNodeBubbles(nodes, elements, boundaries, edgeparas, GRAD);
            //writeMeshRadii(nodes, elements, "meshradii.dat");

            insertNodeBubblesAtBoundaries(boundaries, nodes, elements, edgeparas, GRAD);
            //writeMeshRadii(nodes, elements, "meshradii.dat");

            updatePositionsOfInsertedBoundaryNodeBubbles(nodes, boundaries);
            //writeMeshRadii(nodes, elements, "meshradii.dat");

            updateMeshConnections(edgeparas, elements, nodes);
            //writeMeshRadii(nodes, elements, "meshradii.dat");

            const std::size_t sz_new{ nodes.size() };
            if (sz_new != sz_old)
                round = 0;
            else {
                round++;
                if (round > 2) break;
            }
        }

        //writeMeshRadii(nodes, elements, "meshradii.dat");
        resetRadiiOfBoundaryNodeBubblesAccordingToBoundaryEdgeLength(nodes, elements, edgeparas, GRAD);
        //writeMeshRadii(nodes, elements, "meshradii.dat");

        updateRadiiOfBoundaryNodeBubblesAroundNarrowRegions(nodes, boundaries, edgeparas, GRAD);
        //writeMeshRadii(nodes, elements, "meshradii.dat");

        updateRadiiOfBoundaryNodeBubblesAccordingToGradient(nodes, elements, edgeparas, GRAD);        
        //writeMeshRadii(nodes, elements, "meshradii.dat");

        //Step 2: Mesh Generation for Interior Regions 

        labelBoundaryNodeBubblesAroundNarrowRegions(nodes, boundaries, edgeparas, GRAD);
        //writeMeshRadii(nodes, elements, meshradiifile); 

        generateInteriorMeshes(nodes, elements, edgeparas, GRAD, RMAX);
        writeMeshRadii(nodes, elements, meshradiifile);
    }
};

//int main(){
//    std::ios_base::sync_with_stdio(false);    
//    const std::filesystem::path geofile{ "geometry.dat" };
//    const std::filesystem::path meshfeaturefile{ "meshFeature.txt" };
//    const std::filesystem::path meshradiifile{ "meshradii.dat" };
//    MeshGeneration mg;
//    mg.implement(geofile, meshfeaturefile, meshradiifile);
//}

// meshgeneration "geometry.dat" "meshFeature.txt" "meshradii.dat"
int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    if (argc != 4) {
        std::cerr << "Too mang or few arguments\n";
        return 1;
    }
    const std::filesystem::path geofile{ argv[1] };
    const std::filesystem::path meshfeaturefile{ argv[2] };
    const std::filesystem::path meshradiifile{ argv[3] };
    MeshGeneration mg;
    mg.implement(geofile, meshfeaturefile, meshradiifile);
}