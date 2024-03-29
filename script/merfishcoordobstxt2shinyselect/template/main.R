suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(jlutils))
f='INFILE'
outfile='OUTFILE'
mat=read.txt(f)
head(mat)
str(mat)

ui=shiny::fluidPage(
	shiny::titlePanel("Choose cells for a subset"), 
	shiny::sidebarLayout(
		shiny::sidebarPanel(
			shiny::actionButton("choose_toggle", "Choose/unchoose")
			, shiny::actionButton("reset", "Clear")
			, shiny::actionButton("done", "Done")
			, shiny::h3("Instructions:")
			, shiny::tags$ol(
				shiny::tags$li("Highlight points by clicking and dragging.")
				, shiny::tags$li("Click the 'Choose/unchoose' button.")
				, shiny::tags$li("Repeat until all of the desired cells are black.")
				, shiny::tags$li("Click 'Done'.")
				)
			, shiny::h4("Details:")
			, shiny::tags$ul(
				shiny::tags$li("To start over, click 'Clear'")
				, shiny::tags$li(paste("You can also choose/unchoose specific cells", "by clicking on them directly"))
				)
			)
		, shiny::mainPanel(
			shiny::plotOutput("plot1", height="auto", click="plot1_click", brush=shiny::brushOpts(id="plot1_brush"))
			)
		)
	)
server=function(input, output, session) {
	vals=shiny::reactiveValues(keeprows=rep(FALSE, nrow(mat)))
	output$plot1=shiny::renderPlot({
		plot(mat[, c('center_x', 'center_y')], pch='.')
		points(mat[vals$keeprows, c('center_x', 'center_y')], pch='.', col='red')
	}, height=function() {
		session$clientData$output_plot1_width
	})
	shiny::observeEvent(
		input$plot1_click
		, {
			res=shiny::nearPoints(mat, xvar='center_x', yvar='center_y', input$plot1_click, allRows=TRUE)
			vals$keeprows=vals$keeprows | res$selected_
		})
	shiny::observeEvent(
		input$choose_toggle
		, {
			res=shiny::brushedPoints(mat, xvar='center_x', yvar='center_y', input$plot1_brush, allRows=TRUE)
			vals$keeprows=vals$keeprows | res$selected_
		})
	shiny::observeEvent(input$reset, {
		vals$keeprows=rep(FALSE, nrow(mat))
		})
	shiny::observeEvent(input$done, {
		shiny::stopApp(vals$keeprows)
		})
}
sel=suppressMessages(shiny::runApp(shiny::shinyApp(ui, server)))
print(table(sel))
submat=mat[sel, ]
write.table(submat, file=gzfile(outfile), sep='\t', row.names=F, quote=F, col.names=T)
