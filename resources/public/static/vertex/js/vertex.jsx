var snipPrefix = function(s) {
  return s.substring(s.indexOf(':') + 1);
}

var VertexEdges = React.createClass({
  getInitialState: function() {
    return {};
  },

  render: function() {
    var props = this.props;
    var prefix = props.edges[0].split(':')[0]
    var header = props.label + " (" + props.direction + " " + prefix + ")";

    var items = props.edges.map(gid => (
      <ExpandoItem key={gid}>
        <a onClick={() => props.navigate(gid)}>{snipPrefix(gid)}</a>
      </ExpandoItem>
    ));

    return <Expando header={header}>{items}</Expando>;
  }
});

function PropertyRow(props) {
  return (<tr>
    <td className="prop-key mdl-data-table__cell--non-numeric">{props.name}</td>
    <td className="mdl-data-table__cell--non-numeric">{props.value}</td>
  </tr>)
}

var PropertiesView = function(props) {
  var properties = Object.keys(props.vertex.properties).map(function(key) {
    var v = props.vertex.properties[key];
    return <PropertyRow key={key} name={key} value={v} />
  });

  return (
    <div>
      <div className="vertex-properties">
        <table
          className="prop-table mdl-data-table mdl-js-data-table mdl-data-table--selectable mdl-shad--2dp"
        ><tbody>
          {properties}
        </tbody></table>
      </div>
    </div>
  )
}

var EdgesView = function(props) {
  console.log(props)
  var inEdges = Object.keys(props.vertex['in'])
  // Filter out edges with "hasInstance" in label
  .filter(key => key != 'hasInstance')
  .map(function(key) {
    return <VertexEdges
      key={key}
      label={key}
      navigate={props.navigate}
      edges={props.vertex['in'][key]}
      direction="from"
    />
  });
   var outEdges = Object.keys(props.vertex['out'])
  // Filter out edges with "hasInstance" in label
  .filter(key => key != 'hasInstance')
  .map(function(key) {
    return <VertexEdges
      key={key}
      label={key}
      navigate={props.navigate}
      edges={props.vertex['out'][key]}
      direction="to"
    />
  });

  return (
    <div>
      <div className="vertex-edges-wrapper">
        <div className="vertex-in-edges vertex-edges">
          <h4>In Edges</h4>
          {inEdges}
        </div>
        <div className="vertex-out-edges vertex-edges">
          <h4>Out Edges</h4>
          {outEdges}
        </div>
      </div>
    </div>
  )
}

var VertexInput = React.createClass({
  componentDidMount() {
    componentHandler.upgradeElement(this.refs.mdlWrapper)
  },
  render() {
    return <div
      className="mdl-textfield mdl-js-textfield mdl-textfield--floating-label"
      ref="mdlWrapper"
    >
      <label
        className="mdl-textfield__label"
        htmlFor="vertex-gid-input"
      >Enter a vertex GID</label>
      <input
        id="vertex-gid-input"
        type="text"
        name="gid"
        className="mdl-textfield__input"
        onChange={e => this.props.onChange(e.target.value)}
        value={this.props.value}
      />
    </div>
  },
})


var Expando = React.createClass({
  getInitialState() {
    return {
      collapsed: true,
    }
  },
  componentDidMount() {
    var content = $(this.refs.content)
    content.css('margin-top', -content.height());
  },
  onClick() {
    this.setState({collapsed: !this.state.collapsed})
  },
  render() {
    var props = this.props;
    var rootClassName = classNames("expando", "mdl-collapse", "mdl-navigation", {
      "mdl-collapse--opened": !this.state.collapsed,
    })
    
    return (<div className={rootClassName}>
      <a className="mdl-navigation__link mdl-collapse__button expando-header" onClick={this.onClick}>
        <i className="material-icons mdl-collapse__icon mdl-animation--default">expand_more</i>
        {props.header}
      </a>
      <div className="mdl-collapse__content-wrapper expando-content">
        <div className="mdl-collapse__content mdl-animation--default" ref="content">
          {props.children}
        </div>
      </div>
    </div>)
  }
})

function ExpandoItem(props) {
  return <span className="mdl-navigation__link">{props.children}</span>
}

var PubmedLink = function(props) {
  var url = "https://www.ncbi.nlm.nih.gov/pubmed/" + props.id;
  console.log(url);
  return (<div><a href={url} target="_blank">{url}</a></div>)
}

var VertexViewer = React.createClass({
  getInitialState() {
    return {
      input: this.getGIDFromURL(),
      loading: false,
      error: "",
      vertex: {},
    };
  },

  getGIDFromURL() {
    return getParameterByName("gid")
  },

  componentDidMount() {
    window.onpopstate = this.onPopState
    if (this.state.input) {
      this.setVertex(this.state.input, true)
    }
  },

  onPopState(e) {
    var hash = this.getGIDFromURL();
    if (e.state && e.state.gid) {
      this.setVertex(e.state.gid, true)
    } else if (hash) {
      this.setVertex(hash, true)
    } else {
      this.setVertex()
    }
  },

  setVertex(gid, nopushstate) {
    if (!gid) {
      this.setState({vertex: {}, error: ""})
    } else {
      var url = "/gaia/vertex/find/" + gid;
      this.setState({input: gid, loading: true, error: ""});
      $.ajax({
        url: url,
        dataType: 'json',
        type: 'GET',
        success: result => {
          if (Object.keys(result).length > 0) {
            this.setState({vertex: result, loading: false, error: ""})
            if (!nopushstate) {
              // Only push state to history if we found an actual vertex
              // This avoids pushing state for intermediate queries.
              history.pushState({gid: gid}, "Vertex: " + gid, '?gid=' + gid);
            }
          } else {
            this.setState({vertex: {}, loading: false, error: ""})
          }
        },
        error: (xhr, status, err) => {
          this.setState({loading: false, error: err.toString()})
          console.error(url, status, err.toString())
        },
        timeout: 3000,
      });
    }
  },

  render: function() {
    var loading = "";
    if (this.state.loading) {
      loading = <div className="mdl-spinner mdl-js-spinner is-active"></div>
    }

    var error;
    if (this.state.error) {
      error = <div>Request error: {this.state.error}</div>
    }

    var emptyMessage = "";
    if (this.state.input) {
      emptyMessage = "No vertex found";
    }

    var vertex = <div className="empty-vertex">{emptyMessage}</div>;
    var visualizations = [];

    // The vertex isn't empty, so create a VertexView
    if (this.state.vertex.properties) {
      vertex = (<div><PropertiesView vertex={this.state.vertex} /><EdgesView vertex={this.state.vertex} navigate={this.setVertex} /></div>)

      if (this.state.vertex.properties.type === 'Pubmed') {
        console.log(this.state.vertex.properties);
        var link = (<PubmedLink key="pubmed-link" id={this.state.vertex.properties.pmid} />)
        visualizations.push(link);
      }
    }

    return (
      <div>
        <VertexInput onChange={this.setVertex} value={this.state.input} />
        {loading}
        {visualizations}
        {error}
        {vertex}
      </div>
    );
  }
});

window.addEventListener('load', function() {
  ReactDOM.render(<VertexViewer />, document.getElementById('vertex-explore'));
})

function getParameterByName(name, url) {
    if (!url) {
      url = window.location.href;
    }
    name = name.replace(/[\[\]]/g, "\\$&");
    var regex = new RegExp("[?&]" + name + "(=([^&#]*)|&|#|$)"),
        results = regex.exec(url);
    if (!results) return null;
    if (!results[2]) return '';
    return decodeURIComponent(results[2].replace(/\+/g, " "));
}
