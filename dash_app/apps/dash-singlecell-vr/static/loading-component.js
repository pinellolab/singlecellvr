AFRAME.registerComponent('loading', {
    schema: {
        time: {type: 'number', default: 1000},
        show: {type: 'boolean', default: true},
        zDepth: {type: 'number', default: -1}
    },
    init: function() {
        this.initialized = performance.now();
        this.adjusted = false;
        this.width = this.visibleWidthAtZDepth(this.data.zDepth);
        this.height = this.visibleHeightAtZDepth(this.data.zDepth);
    },
    update: function(oldData) {
        if (this.data.show !== oldData.show && !this.data.show) {
            const dismissTime = performance.now();
            const elapsedSeconds = (dismissTime - this.initialized) / 1000;
            if (elapsedSeconds > this.data.time) {
                this.el.object3D.visible = false;
                document.getElementById("hud").setAttribute('visible', true);
                document.getElementById("player-camera").setAttribute("look-controls", "");
                document.getElementById("scene").setAttribute("vr-mode-ui", "enabled: true");
            } else {
                setTimeout(() => {
                    this.el.object3D.visible = false;
                    document.getElementById("hud").setAttribute('visible', true);
                    document.getElementById("player-camera").setAttribute("look-controls", "");
                    document.getElementById("scene").setAttribute("vr-mode-ui", "enabled: true");
                }, ((this.data.time - elapsedSeconds) * 1000))
            }
        }
        const width = this.visibleWidthAtZDepth(this.data.zDepth);
        const height = this.visibleHeightAtZDepth(this.data.zDepth);
        // if (!this.adjusted) {
        //     console.log("This width: ", this.width);
        //     console.log("This height: ", this.height);
        //     console.log("Width: ", width);
        //     console.log("Height: ", height);
        //     this.width = width;
        //     this.height = height;
        //     this.adjusted = true;
        // }
        const matchMedia = window.webkitMatchMedia || window.mozMatchMedia || window.oMatchMedia || window.msMatchMedia || window.matchMedia;
        const isPortrait = (matchMedia && matchMedia("(orientation: portrait)").matches) || 
                            (["portrait", "portrait-primary", "portrait-secondary"].includes(window.screen.orientation) ||
                            (Utils.mobilecheck() && window.height > window.width))
        if (isPortrait) {
            this.el.object3D.rotation.set(0, 0, -32.987);
            this.el.setAttribute('width', this.height);
            this.el.setAttribute('height', this.width);
            document.getElementById('hud').object3D.position.set(-width/2 + .25, height/2 - .25, this.data.zDepth);
        } else {
            this.el.object3D.rotation.set(0, 0, 0);
            this.el.setAttribute('width', this.height);
            this.el.setAttribute('height', this.width);
            document.getElementById('hud').object3D.position.set(-width/2 + .25, height/2 - .25, this.data.zDepth);
        }
    },
    tick: function() {
        const cameraPosition = document.getElementById('rig').object3D.position;
        this.el.object3D.position.set(cameraPosition.x, cameraPosition.y, cameraPosition.z - 1);
    },
    // These are duplicated, fix later
    visibleHeightAtZDepth: function( depth ) {
        const camera = this.el.sceneEl.camera;
        // compensate for cameras not positioned at z=0
        const cameraOffset = camera.position.z;
        if ( depth < cameraOffset ) depth -= cameraOffset;
        else depth += cameraOffset;
      
        // vertical fov in radians
        const vFOV = camera.fov * Math.PI / 180;
      
        // Math.abs to ensure the result is always positive
        return 2 * Math.tan( vFOV / 2 ) * Math.abs( depth );
    },
    visibleWidthAtZDepth: function( depth ) {
        const camera = this.el.sceneEl.camera;
        const height = this.visibleHeightAtZDepth( depth, camera );
        let width = height * camera.aspect;
        return width;
    }

});