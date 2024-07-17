#=

Drawing Julia using a Fourier series.
A high definition animation can be seen here: https://youtu.be/rrmx2Q3sO1Y

This code is based on code kindly provided by ric-cioffi (https://github.com/ric-cioffi)
But was rewritten for v0.3 by Ole Kröger.
=#

using Javis, FFTW, FFTViews
using TravelingSalesmanHeuristics

function ground(args...)
    Javis.background("white")
    sethue("black")
end

function circ(; r = 20, vec = O, action = :stroke, color = "grey")
    sethue(color)
    circle(O, r, action)
    my_arrow(O, vec)
    return vec
end

function my_arrow(start_pos, end_pos)
    arrow(
        start_pos,
        end_pos;
        linewidth = distance(start_pos, end_pos) / 100,
        arrowheadlength = 7,
    )
    return end_pos
end

function draw_line(
    p1 = O,
    p2 = O;
    color = "grey",
    action = :stroke,
    edge = "solid",
    linewidth = 3,
)
    sethue(color)
    setdash(edge)
    setline(linewidth)
    line(p1, p2, action)
end


using Colors
# , Luxor

# let colors = collect(Colors.color_names)
#     function randomhue()
#         sethue(colors[Int64(ceil(rand() * 100))][1])
#     end
# end

global current_color = "red"

function draw_path!(path, pos, color)
    sethue(color)
    # if (mod(length(path), 71) == 0 || path == 0)
    #     cols = collect(Colors.color_names)
    #     clab = convert(Lab, parse(Colorant, cols[Int64(ceil(rand() * 100))][1]))
    #     labelbrightness = 100 - clab.l
    #     sethue(convert(RGB, Lab(labelbrightness, clab.b, clab.a)))
    #     final_color = convert(RGB, Lab(labelbrightness, clab.b, clab.a))
    #     current_color = cc

    # push!(path, pos)
    # return draw_line.(path[2:end], path[1:(end - 1)]; color = cc)
    # else
    push!(path, pos)
    return draw_line.(path[2:end], path[1:(end - 1)]; color = color)
    # end
end


c2p(c::Complex) = Point(real(c), imag(c))

remap_idx(i::Int) = (-1)^i * floor(Int, i / 2)
remap_inv(n::Int) = 2n * sign(n) - 1 * (n > 0)

function interpol(p1, p2, t)
    return (p2 - p1) * t + p1
end

function param_curve(points, x)
    n = length(points)
    if (x == 1.0)
        return points[n - 1]
    end
    s = trunc(Int, x * (n - 1)) + 1
    e = s + 1
    t = x * (n - 1) - (s - 1)
    return interpol(points[s], points[e], t)
end


using Luxor: polymove!, readsvg
const pagewidth = 1190 # A2 size paper, points
const pageheight = 1684 # A2 size paper, points

function get_points(npoints, options)

    Drawing() # julialogo needs a drawing
    # fontface("DejaVu Sans")
    # textpath("BuddhiLW")
    # textpath("VJ Vacão 🐮")# VJ Vacão 🐄  # Rapadur
    # print(getpath())
    # plist = pathtopoly()
    # julialogo(; action = :path,centered = true)
    #
    logo_path = Luxor.readsvg("./logo-path")
    Javis.pathsvg("logo_path", :path, centered = true)
    logo_path(; action = :path, centered = true)

    shapes = pathtopoly()
    for (n, pgon) in enumerate(shapes)
        randomhue()
        prettypoly(pgon, :stroke, close = true)
        gsave()
        polymove!(pgon, O, Point(-20, 0))
        poly(polysortbyangle(pgon, polycentroid(pgon)), :stroke, close = true)
        grestore()
    end
    new_shapes = shapes[1:end]
    # plot(shapes)
    last_i = 1
    # the circles in the JuliaLogo are part of a single shape
    # this loop creates new shapes for each circle
    # for shape in shapes[1:2]
    #     max_dist = 0.0
    #     for i in 2:length(shape)
    #         d = distance(shape[i - 1], shape[i])
    #         if d > 3
    #             push!(new_shapes, shape[last_i:(i - 1)])
    #             last_i = i
    #         end
    #     end
    # end


    push!(new_shapes, shapes[length(shapes)][last_i:end])
    shapes = new_shapes
    for i in 1:length(shapes)
        shapes[i] .*= options.shape_scale
    end

    total_distance = 0.0
    for shape in shapes
        total_distance += polyperimeter(shape)
    end
    parts = []
    points = Point[]
    start_i = 1
    for shape in shapes
        len = polyperimeter(shape)
        portion = len / total_distance
        nlocalpoints = floor(Int, portion * npoints)
        new_points = [
            Javis.get_polypoint_at(shape, i / (nlocalpoints - 1)) for
            i in 0:(nlocalpoints - 1)
        ]
        append!(points, new_points)
        new_i = start_i + length(new_points) - 1
        push!(parts, start_i:new_i)
        start_i = new_i
    end
    return points, parts
end

using CSV, Tables

function animate_fourier(options)
    npoints = options.npoints
    nplay_frames = options.nplay_frames
    nruns = options.nruns
    nframes = nplay_frames + options.nend_frames

    # obtain points from julialogo
    # points, parts = get_points(npoints, options)

    #
    all_points =
        CSV.File("points1.csv") |>
        Tables.matrix |>
        transpose |>
        x -> Point.(x[1, :], x[2, :])
    # if options.npoints
    # step = Int64(floor(length(all_points) / npoints))
    points = all_points ./ options.shape_scale ## .- Point(300, 100) ## all_points[1:step:end] ./ options.shape_scale
    θ = -(π / 2 - π / 10)
    points =
        map(p -> Point(p.x * cos(θ) - p.y * sin(θ), p.x * sin(θ) + p.y * cos(θ)), points)

    points = (points .- Point(200, -100)) .* -1
    npoints = length(points)
    println("#points: $npoints")
    # solve tsp to reduce length of extra edges
    distmat = [distance(points[i], points[j]) for i in 1:npoints, j in 1:npoints]

    path, cost = solve_tsp(distmat; quality_factor = options.tsp_quality_factor)
    println("TSP cost: $cost")
    points = points[path] # tsp saves the last point again

    # optain the fft result and scale
    x = [p.x for p in points]
    y = [p.y for p in points]

    fs = FFTView(fft(complex.(x, y)))
    # normalize the points as fs isn't normalized
    fs ./= npoints
    npoints = length(fs)

    video = Video(options.width, options.height)
    Javis.Background(1:nframes, ground)

    circles = Object[]

    for i in 1:npoints
        ridx = remap_idx(i)

        push!(circles, Object((args...) -> circ(; r = abs(fs[ridx]), vec = c2p(fs[ridx]))))

        if i > 1
            # translate to the tip of the vector of the previous circle
            act!(circles[i], Action(1:1, anim_translate(O, circles[i - 1])))
        end
        ridx = remap_idx(i)
        act!(circles[i], Action(1:nplay_frames, anim_rotate(0.0, ridx * 2π * nruns)))
    end

    trace_points = Point[]
    Object(1:nframes, (args...) -> draw_path!(trace_points, pos(circles[end]), "red"))

    return render(video; pathname = options.filename)
end

function main()
    hd_options = (
        npoints = 3001, # rough number of points for the shape => number of circles
        nplay_frames = 1200, # number of frames for the animation of fourier
        nruns = 2, # how often it's drawn
        nend_frames = 200,  # number of frames in the end
        width = 1920,
        height = 1080,
        shape_scale = 25, # scale factor for the logo
        tsp_quality_factor = 50,
        filename = "julia_hd.mp4",
    )

    fast_options = (
        npoints = 1001, # rough number of points for the shape => number of circles
        nplay_frames = 600, # number of frames for the animation of fourier
        nruns = 1, # how often it's drawn
        nend_frames = 200,  # number of frames in the end
        width = 1000,
        height = 768,
        shape_scale = 15, # scale factor for the logo
        tsp_quality_factor = 40,
        filename = "julia_fast.mp4",
    )

    gif_options = (
        npoints = 300, # 651, # rough number of points for the shape => number of circles
        nplay_frames = 100, # 600, # number of frames for the animation of fourier
        nruns = 2, # how often it's drawn
        nend_frames = 0,  # number of frames in the end
        width = 1000,
        height = 400,
        shape_scale = 6, # scale factor for the logo
        tsp_quality_factor = 0,
        filename = "julia_logo_dft.gif",
    )

    proportion_options = (
        # # Test
        # npoints = 651, # 651, # rough number of points for the shape => number of circles
        # nplay_frames = 600, # 600, # number of frames for the animation of fourier
        # nruns = 1, # how often it's drawn
        # nend_frames = 3,  # number of frames in the end


        # Render
        npoints = 3000, # rough number of points for the shape => number of circles
        nplay_frames = 1500, # number of frames for the animation of fourier
        nruns = 1, # how often it's drawn
        nend_frames = 30,  # number of frames in the end
        width = 1100,
        height = 700,
        shape_scale = 19, # scale factor for the logo
        tsp_quality_factor = 45,
        filename = "./gifs/prop-blw-logo-color.mp4",
    )
    return animate_fourier(proportion_options)
end

main()
