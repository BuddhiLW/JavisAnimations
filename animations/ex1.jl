using Javis, Animations

# define ground function here

video = Video(500, 500)
translate_anim = Animation(
    [0, 1], # must go from 0 to 1
    [O, Point(150, 0)],
    [sineio()],
)

translate_back_anim = Animation(
    [0, 1], # must go from 0 to 1
    [O, Point(-150, 0)],
    [sineio()],
)

rotate_anim = Animation(
    [0, 1], # must go from 0 to 1
    [0, 2π],
    [linear()],
)

Background(1:150, ground)
ball = Object((args...) -> circle(O, 25, :fill))
act!(ball, Action(1:10, sineio(), scale()))
act!(ball, Action(11:50, translate_anim, translate()))
act!(ball, Action(51:100, rotate_anim, rotate_around(Point(-150, 0))))
act!(ball, Action(101:140, translate_back_anim, translate()))
act!(ball, Action(141:150, rev(sineio()), scale()))

render(video)
