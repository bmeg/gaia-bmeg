var PieChart = React.createClass({
  
  getInitialState: function() {
    return {results:[]}
  },
  
  query: function() {
    Ophion().query().has("gid", ["gene:BRCA1"]).incoming("inGene").groupCount("variantClassification").by("variantClassification").cap(["variantClassification"]).execute(
      function(results) {
        
        var res = []
        for (var k in results['result'][0]) {
          res.push( {"title" : k, "value" : results['result'][0][k]} )
        }
        console.log("res")
        console.log(res)
        this.setState( {"results" : res} )
        this.pieupdate()
      }.bind(this)      
    )
  },
  
  componentDidMount(){
    this.query();
  },

  pieupdate: function() {
    console.log(this.state);
    var cohort = this.state.results;
    
    var pie = d3.layout.pie().value(function(d) {return d.value});
    var slices = pie(cohort);

    var arc = d3.svg.arc().innerRadius(0).outerRadius(100);
    var color = d3.scale.category10();

    var svg = d3.select('svg.piechart');
    var g = svg.append('g').attr('transform', 'translate(200, 200)');

    console.log(svg);
    console.log(slices);

    g.selectAll('path.piechart')
      .data(slices)
      .enter()
      .append('path')
      .attr('class', 'slice')
      .attr('d', arc)
      .attr('fill', function(d) {return color(d.data.value)});

    svg.append('g')
      .attr('class', 'legend')
      .selectAll('text')
      .data(slices)
      .enter()
      .append('text')
      .text(function(d) { return '- ' + d.data.title; })
      .attr('fill', function(d) { return color(d.data.title); })
      .attr('y', function(d, i) { return 20 * (i + 1); })
  },
  
  render: function() {
    return (
      <div></div>
    )
  }  
  
})

ReactDOM.render(<PieChart />, document.getElementById('mutation-types'));
