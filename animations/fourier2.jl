using Javis
using FFTW
using FFTViews
using FileIO
using Images
using TravelingSalesmanHeuristics
using StatsBase

function ground(args...)
    background("black")
    sethue("white")
end

function circ(; r = 10, vec = O, action = :stroke, color = "white")
    sethue("black")
    circle(O, r, action)
    # my_arrow(O, vec)
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
    color = "white",
    action = :stroke,
    edge = "solid",
    linewidth = 3,
)
    sethue(color)
    setdash(edge)
    Luxor.setline(linewidth)
    line(p1, p2, action)
end

function draw_path!(path, pos, color)
    sethue(color)

    push!(path, pos)
    return draw_line.(path[2:end], path[1:(end - 1)]; color = color)
end

function get_points(img)
    # findall(x -> x == 1, img) .|> x -> Point(x.I)
    findall(x -> x !== RGBA{N0f8}(0.0, 0.0, 0.0, 1.0), img) .|> x -> Point(x.I)
end

# function texty()
#     fontsize(36)
#     text("Loading Asset: Jacob Zelko", Point(0, 290); halign = :center)
# end

c2p(c::Complex) = Point(real(c), imag(c))

remap_idx(i::Int) = (-1)^i * floor(Int, i / 2)
remap_inv(n::Int) = 2n * sign(n) - 1 * (n > 0)

using TestImages, CSV, Tables

function animate_fourier(options)
    npoints = options.npoints
    nplay_frames = options.nplay_frames
    nruns = options.nruns
    nframes = nplay_frames + options.nend_frames

    # ImageView.imshow(load("./logo-border2.png"))
    # obtain points from julialogo
    #load(File(format"PNG", "logo-border.png"))

    # all_points = get_points(load("logo-border2.png"))
    # step = Int64(floor(length(all_points) / npoints))
    # points = all_points[1:step:end]

    all_points =
        CSV.File("points0.csv") |>
        Tables.matrix |>
        transpose |>
        x -> Point.(x[1, :], x[2, :])
    step = Int64(floor(length(all_points) / npoints))
    points = all_points[1:step:end] ./ options.shape_scale

    npoints = length(points)
    println("#points: $npoints")
    # println(points)
    # solve tsp to reduce length of extra edges
    distmat = [distance(points[i], points[j]) for i in 1:npoints, j in 1:npoints]

    path, cost = solve_tsp(distmat; quality_factor = options.tsp_quality_factor)
    println("TSP cost: $cost")

    points = points[path] # tsp saves the last point again

    # optain the fft result and scale
    y = [p.x - 372.6666666667 for p in points]
    x = [p.y - 323.3333333333 / 1.25 for p in points]

    fs = FFTView(fft(complex.(x, y)))
    # normalize the points as fs isn't normalized
    fs ./= npoints
    npoints = length(fs)

    video = Video(options.width, options.height)
    Background(1:nframes, ground)

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
    Object(1:nframes, (args...) -> draw_path!(trace_points, pos(circles[end]), "white"))

    # loading = Object(1:nframes, (args...) -> texty())
    # act!(loading, Action(1:400, appear(:draw_text)))

    return render(video; pathname = joinpath(@__DIR__, options.filename))
    # return render(video; liveview = true)
end

function main()
    # hd_options = (
    # npoints = 3001, # rough number of points for the shape => number of circles
    # nplay_frames = 1200, # number of frames for the animation of fourier
    # nruns = 2, # how often it's drawn
    # nend_frames = 200,  # number of frames in the end
    # width = 1920,
    # height = 1080,
    # shape_scale = 2.5, # scale factor for the logo
    # tsp_quality_factor = 50,
    # filename = "julia_hd.mp4",
    # )

    gif_options = (
        npoints = 30, # rough number of points for the shape => number of circles
        nplay_frames = 100, # number of frames for the animation of fourier
        nruns = 1, # how often it's drawn
        nend_frames = 30,  # number of frames in the end
        width = 1100,
        height = 1100,
        shape_scale = 13, # scale factor for the logo
        tsp_quality_factor = 45,
        filename = "./gifs/blw-logo.mp4",
    )

    # gif_options = (
    # npoints = 651, # rough number of points for the shape => number of circles
    # nplay_frames = 600, # number of frames for the animation of fourier
    # nruns = 2, # how often it's drawn
    # nend_frames = 0,  # number of frames in the end
    # width = 350,
    # height = 219,
    # shape_scale = 0.8, # scale factor for the logo
    # tsp_quality_factor = 80,
    # filename = "julia_logo_dft.gif",
    # )
    return animate_fourier(gif_options)
end

main()
