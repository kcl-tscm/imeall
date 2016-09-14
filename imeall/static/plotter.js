// Plotting Routines for Visualizing different potentials.
var plot_lines = function(d,i){
  svg.append("path")
    .attr("class", "line")
    .attr("id", 'tag'+d.key.replace(/\./,''))
    .attr("d", valueline(d.values))
    .style("stroke", function(){return d.color = color(d.key)})
//Draw dots with tooltip information
  svg.selectAll("dot")
    .data(d.values)
	  .enter().append("circle")
	    .attr("class", "dot")
      .attr("id", 'tagdot'+d.key.replace(/\./,''))
	    .attr("r", 3.5)
	    .attr("cx", (d.values, function(d) {return x(d.angle); }))
	    .attr("cy", (d.values, function(d) {return y(d.min_en); }))
      .style("fill", function(){
        return d.color = color(d.key) })
	    .on("mouseover", function(d){
	      tooltip.transition()
	        .duration(200)
	        .style("opacity", 0.9)
	      tooltip.html("<p> or: "+d.or_axis + "</p> <p> bp: "+d.bp+"</p>"
                    +"<p> Angle: " + Math.round(d.angle*100)/100 + "</p>"
                    +"<p> Energy: "+ Math.round(d.min_en*1000)/1000 + " J/m^2 </p>")
	    })
	    .on("mouseout", function(d) {
	      tooltip.transition()
	        .duration(500)
	        .style("opacity",0)
      });

//Add legend/add and remove energetics.
  svg.append("text")
    .attr("x", width)
    .attr("y", margin.top + i*15)
    .attr("class", "legend")
    .style("fill", function(){
      return d.color = color(d.key) })
    .on("click", function(){
      var active = d.active ? false: true;
      newOpacity = active ? 0:1;
      d3.select("#tag"+d.key.replace(/\./,''))
        .transition().duration(100)
        .style("opacity", newOpacity);
      d3.select("#tagdot"+d.key.replace(/\./,''))
        .transition().duration(100)
        .style("opacity", newOpacity);

      d.active = active;
      })
    .text(d.key);
  });
svg.append("g")
  .attr("class", "x axis")
  .attr("transform", "translate(0," + height + ")")
  .call(xAxis);
svg.append("g")
  .attr("class", "y axis")
  .call(yAxis)
  });
</script>
